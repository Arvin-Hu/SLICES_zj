#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/nvidia/lib:/usr/local/nvidia/lib64:$LD_LIBRARY_PATH
# export PATH=/opt/miniconda/bin:/opt/miniconda/envs/umat/bin:$PATH
# export PATH=/opt/conda/bin:/opt/conda/envs/chgnet/bin:$PATH


# 创建slurm配置目录并复制配置文件
mkdir -p /etc/slurm
mkdir -p /etc/slurm-llnl
if [ -f "/crystal/slurm.conf" ]; then
    cp /crystal/slurm.conf /etc/slurm/
    cp /crystal/slurm.conf /etc/slurm-llnl/
    chown slurm:slurm /etc/slurm/slurm.conf /etc/slurm-llnl/slurm.conf
    ls /etc/slurm-llnl
    ls /etc/slurm
fi

# echo "127.0.0.1 workq" >> /etc/hosts
# hostname

service munge restart
service slurmctld restart
service slurmd restart
sinfo


source /opt/miniconda/etc/profile.d/conda.sh
conda activate umat #chgnet
cd MatterGPT
echo "Starting Gradio app..."
python app.py
