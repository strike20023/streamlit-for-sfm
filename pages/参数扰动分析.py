import streamlit as st
import time
from pathlib import Path
from utils.run_command import experimental_results_img
from utils.make_plot import display_plot
from 定位测姿误差分析平台 import 加载图片, render_svg, check_button_clicked
import os
import json
st.markdown(
    """<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>""",
    unsafe_allow_html=True,
)

def 加载结果图片(paths, col_n=3):
    # paths = list(Path(paths).rglob("*"))
    # cols = [st.columns(col_n) for _ in range(len(paths) // col_n + 1)]
    # for i, p in enumerate(paths):
    #     cols[i // col_n][i % col_n].image(
    #         p.__str__(), use_column_width=True, caption=p.stem
    #     )

    paths = list(Path(paths).rglob("*.svg"))
    for i, p in enumerate(paths):
        print(p)
        render_svg(p)


def analysis(参数):
    参数['output_excel_path'] = '/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_img/excel/test.xlsx'
    参数['output_img_dir'] = '/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_img/img'
    # os.system(f'rm -rf {param['output_img_dir']}/*')

    with st.status("正在分析", expanded=False) as status:
        status.update(label="正在读取图片...")
        experimental_results_img(参数, status, st)
        # status.update(label="正在重建结果...")
        加载结果图片(参数['output_img_dir'])
        # status.update(label="分析结束", state="complete", expanded=True)


def app():
    加载图片("data")
    st.write("#### 参数设置")
    param = dict()
    图片组 = st.selectbox("选择图片组", [i.stem for i in list(Path("data").glob('*')) if (i/'images').exists()])
    col = st.columns(4)
    to_str = lambda x: ' '.join([str(i) for i in x])
    param.update(
        {
            'input_img': 图片组,
            'variance': to_str(col[1].slider("特征点偏移方差范围", 0.1, 3.0, (0.3, 2.1), 0.05)),
            'percentage': to_str(col[0].slider("特征点偏移百分比范围", 0., 1.0, (0., 1.), 0.05)),
            'fxy': to_str(col[2].slider("焦距偏移范围", -.6, .6, (-0.3, 0.3), 0.05)),
            'cxy': to_str(col[3].slider("光轴偏移范围", -.6, .6, (-0.3, 0.3), 0.05)),
        }
    )

    print(param)
    if check_button_clicked('开始计算'):
        analysis(param)
    if param['variance'] == '0.3 2.1' \
        and param['percentage'] == '0.0 1.0' \
        and param['fxy'] == '-0.3 0.3' \
        and param['cxy'] == '-0.3 0.3':
        # 分析_btn = st.button("开始分析")
        
        if check_button_clicked("开始分析"):
            # 保存状态
            display_plot(st, 'result_img/eg')

app()
