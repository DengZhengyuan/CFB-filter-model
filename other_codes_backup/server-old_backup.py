# -*- coding: utf-8 -*-
print("Let's start!")
# ---------------------------------------------------------
# socket_echo_server_uds.py
import sys
import socket
import os
import struct
import time
import traceback
import threading
import multiprocessing
from multiprocessing import Process


os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
print("Import packages and disable GPU complete. ")
# ---------------------------------------------------------
import numpy as np
from sklearn.preprocessing import StandardScaler
import pickle

print("Import pickle complete. ")
# ----------------------------------

model_hr_address = "./model/model_3D-hr-cubic-half.h5"
model_hd_address = "./model/model_3D-hd-cubic-half.h5"

# ----------------------------------
# for Hr
# 0 nt
# 1 avEs
# 2 gradXEs
# 3 gradYEs
# 4 gradXP
# 5 gradYP
# 6 avUslip
# 7 gradXUslip
# 8 gradYUslip
# 9 avOz
# 10 gradXOz
# 11 gradYOz
# 12 coefHr
# ----------------------------------

nml_factor_hr = np.array([[
    3.6799843384E+01,
    1.5982814355E-01,
    1.2721378852E+03,
    4.2235212576E-01,
    1.1953221242E+06,
    8.8954894382E+00,
    3.2707404349E-01,
    8.7746568175E+03,
    2.6247173478E-01,
    5.6427935728E-01,
    2.2325880651E+03,
    1.4488878672E+00,
    6.3149212593E-01
],
[ 
    1.1704891037E+01,
    1.1823957392E-01,
    2.2825528724E+03,
    4.0476079289E+01,
    2.8380816442E+06,
    1.9503716592E+03,
    2.4506833706E-01,
    1.9736931154E+04,
    7.6984981100E+01,
    2.8701911491E-01,
    5.1183006909E+03,
    5.5640726584E+01,
    7.2881395757E-02
]])

# ----------------------------------
# for Hd
# 0 nt
# 1 avEs
# 2 gradXEs
# 3 gradYEs
# 4 gradXP
# 5 gradYP
# 6 avUslipY
# 7 gradXUslipY
# 8 gradYUslipY
# 9 avBetay
# 10 coefHd
# ----------------------------------

nml_factor_hd = np.array([[
    3.70035695213E+01,
    1.58755329696E-01,
    1.29003760913E+03,
    9.75815987820E-01,
    1.15935057931E+06,
    2.73008873578E+01,
    3.09356620558E-01,
    8.13747803847E+03,
    2.92646374665E-01,
    3.34050820274E+04,
    5.21603190897E-01
],
[
    1.18840855360E+01,
    1.17953784620E-01,
    2.31289437630E+03,
    3.98219136394E+01,
    2.77709590070E+06,
    1.89276143582E+03,
    2.33254151638E-01,
    1.82825139016E+04,
    6.93896396324E+01,
    4.37588727385E+04,
    2.22624166564E-01
]
])
# ------------------------------
# ---------------------------------------------------------
from tensorflow import keras as keras
from tensorflow.keras.models import load_model
import tensorflow as tf
print("Import tensorflow complete. ")
# ---------------------------------------------------------
# Set the server address
# The address of socket, same as the address in C codes
server_address = "/tmp/welfkewgsocket"
# ---------------------------------------------------------
# Make sure the socket does not already exist
# 判断是否socket已经被占用
# 如果try，则执行except
try:
    os.unlink(server_address)
except OSError:
    print("OSError")
    if os.path.exists(server_address):
        raise

# Create a UDS socket
# 创建socket， 不用改
sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)

# Bind the socket to the address
print("starting up on {}".format(server_address))
sock.bind(server_address)  # socket 的地址

# Listen for incoming connections
# socket读取数据，有数据才继续，没数据就等着
sock.listen(1)

PACKETSEP = b"\0"  # 分隔符，与c一致
KILL = 0  # kill server # commd的数据，除了数值还会带上
FEATURE = 0  # 特征数据
RESULT = 0  # 返回值

# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
def writePacket(connection, command, buf):  # 不用改
    bufLen = len(buf)
    # Write Header
    try:
        connection.send(PACKETSEP)
        connection.send(struct.pack("b", command))
        connection.send(struct.pack("i", bufLen))
        connection.send(buf)
    except OSError as e:
        print(e)
        print("OSError, client may have disconnected")
        connection.close()
        return -1
    except BrokenPipeError as e:
        print(e)
        print("BrokenPipeError, client may have disconnected")
        connection.close()
        return -1
    except ConnectionResetError as e:
        print(e)
        print("Connection reset by peer")
        connection.close()
        return -1
    return bufLen


# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
def readPacket(connection):
    try:
        # 1. Read till a PACKETSEP
        while True:
            curByte = connection.recv(1)
            if curByte == b"":
                return None, None, None
            if curByte == PACKETSEP:
                break
            pass
        # 2. Read command
        cmd = struct.unpack("b", connection.recv(1, socket.MSG_WAITALL))[0]

        # 3. Read packet length
        packetLen = struct.unpack("i", connection.recv(4, socket.MSG_WAITALL))[0]
        # 4. Read Data
        data = connection.recv(packetLen, socket.MSG_WAITALL)
    except OSError as e:
        print(e)
        print("OSError, client may have disconnected")
        connection.close()
        return None, None, None
    except BrokenPipeError as e:
        print(e)
        print("BrokenPipeError, client may have disconnected")
        connection.close()
        return None, None, None

    return cmd, packetLen, data  # 返回的命令，数据长度，数据


# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
def readPacket1(connection, structStr, errorMsg):  # structstr 需要学习如何写
    cmd, packetLen, data = readPacket(connection)
    # print('readPacket1',cmd,packetLen,data)
    if cmd is None or packetLen is None or data is None:
        return None
    # print('packetLen==struct.calcSize?',packetLen == struct.calcsize(structStr))
    if packetLen == struct.calcsize(structStr):
        # print("Unpacking tuple!!", packetLen, len(data))
        unPackedTuple = struct.unpack(structStr, data)  # 判断数据格式，填进指定的数
        return unPackedTuple
    else:
        print(errorMsg)
        return None


# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------
def mp_worker(connection, client_address, workerId):  # 多进程执行
    print("Worker Started", workerId)
    """
    with open('./ML_Models/standardscaler20190728.bin','rb') as f:           #打开scaler的参数, 归一化参数
    scaler=pickle.load(f)
    """
    # 加载model
    regressorhd = tf.keras.models.load_model(model_hd_address, compile=False)
    regressorhr = tf.keras.models.load_model(model_hr_address, compile=False)

    while True:  # 死循环，worker不断读取c的数据，写入
        # 调用readpacket1，数据长度放入arrLen，‘i’需要按实际类型修改
        arrLen = readPacket1(connection, "i", "Unsupported Data, expect array length")
        if arrLen is None:
            break
        arrLen = arrLen[0]  # 数组中第一个元素

        if arrLen == 0:
            print("Array len is zero, close this process")
            break

        # connection标识符，arrlen：特征个数
        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        arr_nt = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect nt Array"
        )
        if arr_nt is None:
            break
        arr_nt_hr = (np.array(arr_nt) - nml_factor_hr[0,0]) / nml_factor_hr[1,0]
        arr_nt_hd = (np.array(arr_nt) - nml_factor_hd[0,0]) / nml_factor_hd[1,0]
        # print("arr_nt")
        # ------------------------------------------------------------------------------
        arr_Es = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Es Array"
        )
        if arr_Es is None:
            break
        arr_Es_hr = (np.array(arr_Es) - nml_factor_hr[0,1]) / nml_factor_hr[1,1]
        arr_Es_hd = (np.array(arr_Es) - nml_factor_hd[0,1]) / nml_factor_hd[1,1]
        # print("arr_Es")
        # ------------------------------------------------------------------------------
        arr_Es_dx = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Es_dx Array"
        )
        if arr_Es_dx is None:
            break
        arr_Es_dx_hr = (np.array(arr_Es_dx) - nml_factor_hr[0,2]) / nml_factor_hr[1,2]
        arr_Es_dx_hd = (np.array(arr_Es_dx) - nml_factor_hd[0,2]) / nml_factor_hd[1,2]
        # print("arr_Es_dx")
        # ------------------------------------------------------------------------------
        arr_Es_dy = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Es_dy Array"
        )
        if arr_Es_dy is None:
            break
        arr_Es_dy_hr = (np.array(arr_Es_dy) - nml_factor_hr[0,3]) / nml_factor_hr[1,3]
        arr_Es_dy_hd = (np.array(arr_Es_dy) - nml_factor_hd[0,3]) / nml_factor_hd[1,3]
        # print("arr_Es_dy")
        # ------------------------------------------------------------------------------
        arr_P_dx = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect p_dx Array"
        )
        if arr_P_dx is None:
            break
        arr_P_dx_hr = (np.array(arr_P_dx) - nml_factor_hr[0,4]) / nml_factor_hr[1,4]
        arr_P_dx_hd = (np.array(arr_P_dx) - nml_factor_hd[0,4]) / nml_factor_hd[1,4]
        # print("arr_P_dx")
        # ------------------------------------------------------------------------------
        arr_P_dy = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect p_dy Array"
        )
        if arr_P_dy is None:
            break
        arr_P_dy_hr = (np.array(arr_P_dy) - nml_factor_hr[0,5]) / nml_factor_hr[1,5]
        arr_P_dy_hd = (np.array(arr_P_dy) - nml_factor_hd[0,5]) / nml_factor_hd[1,5]
        # print("arr_P_dy")
        # ------------------------------------------------------------------------------
        arr_Uslip = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Uslip Array"
        )
        if arr_Uslip is None:
            break
        arr_Uslip_hr = (np.array(arr_Uslip) - nml_factor_hr[0,6]) / nml_factor_hr[1,6]
        arr_Uslip_hd = (np.array(arr_Uslip) - nml_factor_hd[0,6]) / nml_factor_hd[1,6]
        # print("arr_Uslip")
        # ------------------------------------------------------------------------------
        arr_Uslip_dx = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Uslip_dx Array"
        )
        if arr_Uslip_dx is None:
            break
        arr_Uslip_dx_hr = (np.array(arr_Uslip_dx) - nml_factor_hr[0,7]) / nml_factor_hr[1,7]
        arr_Uslip_dx_hd = (np.array(arr_Uslip_dx) - nml_factor_hd[0,7]) / nml_factor_hd[1,7]
        # print("arr_Uslip_dx")
        # ------------------------------------------------------------------------------
        arr_Uslip_dy = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Uslip_dy Array"
        )
        if arr_Uslip_dy is None:
            break
        arr_Uslip_dy_hr = (np.array(arr_Uslip_dy) - nml_factor_hr[0,8]) / nml_factor_hr[1,8]
        arr_Uslip_dy_hd = (np.array(arr_Uslip_dy) - nml_factor_hd[0,8]) / nml_factor_hd[1,8]
        # print("arr_Uslip_dy")
        # ------------------------------------------------------------------------------
        arr_Oz = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Oz Array"
        )
        if arr_Oz is None:
            break
        arr_Oz = (np.array(arr_Oz) - nml_factor_hr[0,9]) / nml_factor_hr[1,9]
        # print("arr_Oz")
        # ------------------------------------------------------------------------------
        arr_Oz_dx = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Oz_dx Array"
        )
        if arr_Oz_dx is None:
            break
        arr_Oz_dx = (np.array(arr_Oz_dx) - nml_factor_hr[0,10]) / nml_factor_hr[1,10]
        # print("arr_Oz_dx")
        # ------------------------------------------------------------------------------
        arr_Oz_dy = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect Oz_dy Array"
        )
        if arr_Oz_dy is None:
            break
        arr_Oz_dy = (np.array(arr_Oz_dy) - nml_factor_hr[0,11]) / nml_factor_hr[1,11]
        # print("arr_Oz_dy")
        # ------------------------------------------------------------------------------
        arr_beta = readPacket1(
            connection, "%sd" % arrLen, "Unsupported Data, expect beta Array"
        )
        if arr_beta is None:
            break
        arr_beta = (np.array(arr_beta) - nml_factor_hd[0,9]) / nml_factor_hd[1,9]
        # print("arr_beta")
        print("Get all variables.")
        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        # Make prediction
        curSamplehd = np.concatenate(
            [
                arr_nt_hd.reshape(-1, 1),
                arr_Es_hd.reshape(-1, 1),
                arr_Es_dx_hd.reshape(-1, 1),
                arr_Es_dy_hd.reshape(-1, 1),
                arr_P_dx_hd.reshape(-1, 1),
                arr_P_dy_hd.reshape(-1, 1),
                arr_Uslip_hd.reshape(-1, 1),
                arr_Uslip_dx_hd.reshape(-1, 1),
                arr_Uslip_dy_hd.reshape(-1, 1),
                arr_beta.reshape(-1, 1)
            ],
            axis=1
        )

        curSamplehr = np.concatenate(
            [
                arr_nt_hr.reshape(-1, 1),
                arr_Es_hr.reshape(-1, 1),
                arr_Es_dx_hr.reshape(-1, 1),
                arr_Es_dy_hr.reshape(-1, 1),
                arr_P_dx_hr.reshape(-1, 1),
                arr_P_dy_hr.reshape(-1, 1),
                arr_Uslip_hr.reshape(-1, 1),
                arr_Uslip_dx_hr.reshape(-1, 1),
                arr_Uslip_dy_hr.reshape(-1, 1),
                arr_Oz.reshape(-1, 1),
                arr_Oz_dx.reshape(-1, 1),
                arr_Oz_dy.reshape(-1, 1)
            ],
            axis=1
        )

        hd = regressorhd.predict(curSamplehd)
        hr = regressorhr.predict(curSamplehr)

        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------
        hdBytes = struct.pack(str(arrLen) + "d", *hd)  # 转化为c的数据类型
        wrlt = writePacket(connection, 8, hdBytes)  # 写作packet

        hdBytes2 = struct.pack(str(arrLen) + "d", *hr)  # 转化为c的数据类型
        wrlt = writePacket(connection, 8, hdBytes2)  # 写作packet

        if wrlt < 0:
            print("Write Failed")
            break
            print("Worker Closed", workerId)
    connection.close()


workerId = 0
while True:
    # Wait for a connection
    print("Waiting for a connection")
    connection, client_address = sock.accept()
    print("New connection")
    p = Process(target=mp_worker, args=(connection, client_address, workerId))
    p.start()
    workerId += 1
