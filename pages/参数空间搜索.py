import streamlit as st
import time
from pathlib import Path
import base64
import re
import json
import os
from utils.tree_static import get_tree

from 定位测姿误差分析平台 import 加载图片, render_svg
from utils.run_command import RL_inference

def analysis(param):
    param['results_path']=Path('result_tree/eg/results_6-5.csv')
    with st.status("正在分析", expanded=False) as status:
        status.update(label="正在读取图片...")
        # param = RL_inference(param, status)
        status.update(label="正在重建结果...")
        get_tree(
            'utils/分支.csv', 
            param['results_path'],
            'svg'
            )
        status.update(label="分析结果👇", state="complete", expanded=True)
        render_svg(param['results_path'].parent/'tree.svg')

def app():
    加载图片("data")
    st.write("## 参数设置")
    param = dict()
    param['DATA_DIR'] = st.selectbox("选择图片组", sorted([i.stem for i in list(Path("data").glob('*')) if (i/'images').exists()]))
    param.update(
        {
            'reward_threshold':st.number_input("位姿点偏移误差", 0., 1., 0.05, 0.01),
        }
    )
    col = st.columns(4)
    param.update(
        {
            **{ ('percentage_'+['min','max'][idx]):i for idx, i in enumerate(col[0].slider("特征点偏移百分比范围", 0., 1.0, (0., 1.), 0.05))},
            **{ ('variance_'+['min','max'][idx]):i for idx, i in enumerate(col[1].slider("特征点偏移方差范围", 0.1, 3.0, (0.3, 2.1), 0.05))},
            **{ ('fxy_action_'+['min','max'][idx]):i for idx, i in enumerate(col[2].slider("焦距偏移范围", -.6, .6, (-0.6, 0.6), 0.05))},
            **{ ('cxy_action_'+['min','max'][idx]):i for idx, i in enumerate(col[3].slider("光轴偏移范围", -.6, .6, (-0.6, 0.6), 0.05))},
        }
    )
    col = st.columns(4)
    param.update(
        {
            'variance':col[1].number_input("特征点偏移方差步长", 0., 1., 0.1, 0.01),
            'percentage':col[0].number_input("特征点偏移百分比步长", 0., 1., 0.05, 0.01),
            'fxy_action':col[2].number_input("焦距偏移步长", 0., 1., 0.02, 0.01),
            'cxy_action':col[3].number_input("光轴偏移步长", 0., 1., 0.05, 0.01),
        }
    )
    计算_btn = st.button("开始计算")
    if 计算_btn:
            analysis(param)
    info = {i.parent.stem:json.load(open(i, 'r')) for i in Path('result_tree/exp').rglob('info.json')}
    print(info)
    print(param)

    exist_path = None
    for k, v in info.items():
        is_same = True
        for key in param.keys():
            if param[key] != v[key]:
                is_same = False
                break
        if is_same:
            exist_path = k
            break
    if exist_path:
        分析_btn = st.button("开始分析")
        if 分析_btn:
            svg_path = Path('result_tree/exp')/exist_path/'tree.svg'
            if not os.path.exists(svg_path):
                get_tree(
                    'utils/分支.csv', 
                    list((Path('result_tree/exp')/exist_path).glob('*.csv'))[0],
                    'svg'
                    )
            render_svg(svg_path)


app()
