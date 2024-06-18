import streamlit as st
import base64
import os
import base64
from pathlib import Path
import re
st.set_page_config(page_title="定位测姿误差分析平台", layout='wide')

def check_button_clicked(button_name):
    if st.button(button_name):
        st.session_state[button_name] = True
    if button_name not in st.session_state:
        st.session_state[button_name] = False
    return st.session_state[button_name]

@st.cache_resource()
def 加载图片(paths):
    with st.expander("照片组👇",expanded=False):
        for 图片文件夹 in list(Path(paths).glob('*')):
            if not (图片文件夹/'images').exists():
                continue
            图片列表 = list((图片文件夹/'images').glob('*.jpg'))
            cols = st.columns(len(图片列表)+2)
            cols[1].write(图片文件夹.stem)
            for i, p in enumerate(图片列表):
                cols[i+2].image(p.__str__(), use_column_width=True)
            break

# @st.cache_resource()
def render_svg(svg_path):
    with open(svg_path, "r") as f:
        svg = f.read()
    # width="484pt" height="116pt"
    rule = re.compile(r'width="([0-9]+)pt"')
    svg = rule.sub('width="100%"', svg)
    rule = re.compile(r'height="([0-9]+)pt"')
    svg = rule.sub('height="auto"', svg)
    b64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

# def render_svg(svg_path):
#     # with open(svg_path, "r") as f:
#     #     svg = f.read()
#     st.image(svg_path, format='svg', use_column_width=True)


def get_binary_file_downloader_html(bin_file):
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">📝</a>'
    return href

# st.markdown("""<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>""", unsafe_allow_html=True)

st.write("# 欢迎使用定位试姿误差分析平台！")

st.sidebar.success("在上方选择需要进行的实验☝️")

st.markdown(
    """
    定位试姿误差分析平台是一个用于分析定位测试姿误差的平台。

    **👈 从侧边栏选择实验类型**
"""
)
table = [
    [[
    "0.1","0.1","0.1", get_binary_file_downloader_html("pages/参数空间搜索.py"),
]],
    [[
    "0.1","0.1","0.1","0.1", get_binary_file_downloader_html("pages/参数空间搜索.py"),
]],

]

param=[['特征点偏移量',
'焦距偏移量',
'光轴偏移量'],
[
'特征点偏移百分比范围',
'特征点偏移方差范围',
'焦距偏移范围',
'光轴偏移范围'
]]
# # table 表格
# st.write("## 历史分析结果")
# col = st.columns(2)
# col[0].markdown(
#     '参数扰动分析结果\n\n| '+'|'.join([i for i in param[0]])+' | 下载 |\n|'+'|'.join(['-------']*4)+'|\n'+f"\n".join([f"| {'|'.join([t for t in i])} | " for i in table[0]]),
#     unsafe_allow_html=True,
# )
# col[1].markdown(
#     '参数空间搜索结果\n\n| '+'|'.join([i for i in param[1]])+' | 下载 |\n|'+'|'.join(['-------']*5)+'|\n'+f"\n".join([f"| {'|'.join([t for t in i])} | " for i in table[1]]),
#     unsafe_allow_html=True,
# )
# plot
import pandas as pd
df = pd.DataFrame({
  '特征点偏移量': [1, 2, 3],
  '焦距偏移量': [1, 2, 3],
  '光轴偏移量': [1, 2, 3],
  '特征点偏移百分比范围': [1, 2, 3]
})
st.scatter_chart(data=df, x='特征点偏移量', y='焦距偏移量', color='光轴偏移量', size='特征点偏移百分比范围', width=0, height=0, use_container_width=True)
