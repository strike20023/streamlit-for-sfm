import subprocess as sp
import streamlit as st
import re

def RL_inference(args, st_status):
    cmd = ['sudo', '/mnt/HUS728/miniconda3/envs/zhenhanbai/bin/python', '/mnt/ST8000/zhenhanbai/RL/main.py']
    args['DATA_DIR'] = '/mnt/ST8000/huangzhe/海纹项目/streamlit_app/data/{}'.format(args['DATA_DIR'])
    
    cmd += [t for k,v in args.items() for t in ['--'+k,str(v)]]
    st_status.update(label="正在执行强化学习...")
    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True) as proc:
        # 实时输出标准输出
        for line in proc.stdout:
            st_status.update(label=line.strip())
        
        # 实时输出标准错误
        for line in proc.stderr:
            st_status.update(label=line.strip(), state="error")

def experimental_results_img(args, st_status, st):
    cmd = ['sudo', '/mnt/HUS728/miniconda3/envs/py36/bin/python', '/mnt/ST8000/jialei/SFM_E/API/run_sfm.py']
    cmd += [t for k,v in args.items() for t in ['--'+k,str(v)]]
    st_status.update(label="正在执行分析...")
    print(' '.join(cmd))
    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True) as proc:
        # 实时输出标准输出
        for line in proc.stdout:
            rule = re.compile(r'{(.*?)%}')
            if rule.findall(line):
                st_status.update(label=rule.findall(line)[0])

        # 实时输出标准错误
        for line in proc.stderr:
            rule = re.compile(r'{(.*?)%}')
            if rule.findall(line):
                st_status.update(label=rule.findall(line)[0], state="error")
                # 进度条
                st.progress(int(rule.findall(line)[0]))

