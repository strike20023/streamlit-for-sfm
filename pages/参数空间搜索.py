import streamlit as st
import time
from pathlib import Path
import base64
import re
from 定位测姿误差分析平台 import 加载图片, render_svg
from utils.run_command import RL_inference
from utils.make_treev2 import 绘制决策树

st.markdown(
    """<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>""",
    unsafe_allow_html=True,
)



def analysis(param):
    param['results_path']='/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_trees/rh_result.csv'
    with st.status("正在分析", expanded=False) as status:
        status.update(label="正在读取图片...")
        param = RL_inference(param, status)
        status.update(label="正在重建结果...")
        绘制决策树(param['results_path'],param['results_path'].replace('.csv','.dot'))
        status.update(label="分析结果👇", state="complete", expanded=True)
        render_svg(param['results_path'].replace('.csv','.dot')+'.svg')


def app():
    加载图片("data")
    st.write("## 参数设置")
    param = dict()
    param['DATA_DIR'] = st.selectbox("选择图片组", [i.stem for i in list(Path("data").glob('*')) if (i/'images').exists()])
    param.update(
        {
            'reward_threshold':st.number_input("位姿点偏移误差", 0., 1., 0.05, 0.01),
        }
    )
    col = st.columns(4)
    param.update(
        {
            **{ ('percentage_'+['min','max'][idx]):i for idx, i in enumerate(col[0].slider("特征点偏移百分比范围", 0.1, 3.0, (0.3, 2.1), 0.05))},
            **{ ('variance_'+['min','max'][idx]):i for idx, i in enumerate(col[1].slider("特征点偏移方差范围", 0., 1.0, (0., 1.), 0.05))},
            **{ ('fxy_action_'+['min','max'][idx]):i for idx, i in enumerate(col[2].slider("焦距偏移范围", -.6, .6, (-0.3, 0.3), 0.05))},
            **{ ('cxy_action_'+['min','max'][idx]):i for idx, i in enumerate(col[3].slider("光轴偏移范围", -.6, .6, (-0.3, 0.3), 0.05))},
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
    # import json
    # print(json.dumps(param))
    计算_btn = st.button("开始计算")
    if 计算_btn:
            analysis(param)
    if param['variance'] == 0.1 \
        and param['percentage'] == 0.05 \
        and param['fxy_action'] == 0.02 \
        and param['cxy_action'] == 0.05:
        分析_btn = st.button("开始分析")
        if 分析_btn:
            render_svg('/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_tree/eg/四因素_0319.dot.svg')


app()
