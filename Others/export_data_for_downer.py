# Description: This script is used to test the pyfluent package
import ansys.fluent.core as pyfluent
import os
import pandas as pd
import numpy as np
import glob

# Change the working directory
root_path = r"E:\\OneDrive - The University of Western Ontario\\CFD_cases_OD\\2D_downer"
work_path = r"\\W5200-5SIn\\Wave 4 - fric"
work_path = root_path + work_path
os.chdir(work_path)
print("The directory is", os.getcwd())


# Define the case file name
case_file_name = "2D-D-W5200-5GIn-28"

# Define the export name: 1~time~2
export_name_1 = [
    ".\\data_trans_5SIn\\[data trans overall ",
    ".\\data_trans_5SIn\\[data trans radial ",
    ".\\data_trans_5SIn\\[data trans axial ",
]
# ~time~
export_name_2 = "s] 2D-D-W5200-5SIn.csv"

# Define the list of variables and surface IDs
list_fields = [
    "x-coordinate",
    "y-coordinate",
    "u_slip",
    "rr",
    "oz",
    "gas-y-velocity",
    "solids-y-velocity",
    "solids-vof",
]
list_surface_ids = [
    [*range(10, 73, 1)],  # overall
    [*range(73, 82, 1)],  # radial
    [*range(82, 93, 1)],  # axial
]


# Define the initial time and interval time of reading Fluent data files
init_time = 12
interval_time = 1

# Launch Fluent in solver mode
solver_session = pyfluent.launch_fluent(
    mode="solver", precision="double", version="2d", processor_count=8 
)


# --------------------------------------------------------------------------------------------------
# Define Functions
# --------------------------------------------------------------------------------------------------
def get_data_in_one_datafile(field_data, list_fields: "list", list_surface_ids: "list"):
    """
    This function is used to get the data in one datafile
    Input:
        field_data: the field data from the solver session
        list_fields: the list of variables, such as temperature, pressure, etc.
                     They can be checked by function
                     `field_data.get_scalar_field_data.field_name.allowed_values()`
        list_surface_ids: the list of surface IDs
    """
    transaction = solver_session.field_data.new_transaction()
    for name_field in list_fields:
        transaction.add_scalar_fields_request(
            surface_ids=list_surface_ids,
            field_name=name_field,
            node_value=False,
            boundary_value=True,
        )

    payload_data = transaction.get_fields()
    export_data = payload_data[
        (("type", "scalar-field"), ("dataLocation", 1), ("boundaryValues", True))
    ]

    for idx, items in enumerate([*export_data.keys()]):
        if idx == 0:
            df_field_data = pd.DataFrame(export_data[items])
        else:
            df_field_data = pd.concat([df_field_data, pd.DataFrame(export_data[items])])

    return df_field_data


def export_csv(field_data, time, export_name_1, export_name_2, list_surface_ids):
    """
    This function is used to export the data to csv file
    Input:
        field_data: the field data from the solver session
        time: the time of the data
        export_name_1: the export name before the time
        export_name_2: the export name after the time
        list_surface_ids: the list of surface IDs
    """
    # get the data matrix in one datafile
    temp_writing_data = get_data_in_one_datafile(
        field_data, list_fields, list_surface_ids
    )
    # write the data to csv file
    csv_name = export_name_1 + str(time) + export_name_2
    temp_writing_data.to_csv(csv_name, mode="w", header=list_fields, index=False)


# --------------------------------------------------------------------------------------------------
# Main Function
# --------------------------------------------------------------------------------------------------
# Check the health of the Fluent solver
solver_session.health_check_service.is_serving
# Get the TUI session
tui = solver_session.tui

# Read the case file
solver_session.file.read(file_type="case", file_name=case_file_name)

dir_datafile = glob.glob(os.path.join(".\\", "*0.dat.h5"))
print("The data files are:")
for filename in dir_datafile:
    print(filename)

for idx_filename, file_dir in enumerate(dir_datafile):  # iterate over the data files
    # get the time from the file name
    time = float(file_dir.split("-")[-1].partition(".dat")[0])
    if time >= init_time and idx_filename % interval_time == 0:
        filename = file_dir.partition(".dat")[0].partition(".\\")[2]
        # read the data file
        solver_session.file.read(file_type="data", file_name=filename)
        # get the field data
        field_data = solver_session.field_data
        for i in range(0, 1, 1):  # iterate over the export names and surface IDs
            export_csv(
                field_data, time, export_name_1[i], export_name_2, list_surface_ids[i]
            )
            print("Writing the", i, "(th)", "data at time = ", time, "s is done!")

# Exit the solver session (Fluent)
solver_session.tui.exit
print("The solver session is closed!")
