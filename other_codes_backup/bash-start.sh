salloc --account=def-zhangc --cpus-per-task=8 --mem=32G --time=1:00:00
source ~/tensorflow/bin/activate
python -u ./server.py &> output-py.txt &

fluent 3ddp -t 4 -mpi=intel -affinity=0 -ssh -g 
fluent 3ddp -t 8 -mpi=intel -affinity=0 -ssh -g
fluent 3ddp -t 8 -mpi=intel -affinity=0 -ssh -g -i jou-fluent-start.jou