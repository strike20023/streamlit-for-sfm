import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

st.markdown("""<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>""", unsafe_allow_html=True)

@st.cache_data
def scatter_chart(df, title=''):
    fig = px.scatter(df)
    fig.add_hline(y=0.05, line_dash="dash", line_color="red")
    # title
    fig.update_layout(
        title={
            'text': title,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title="",
        yaxis_title="误差百分比",
        # font=dict(
        #     family="Courier New, monospace",
        #     size=18,
        #     color="RebeccaPurple"
        # )
    )
    return fig
    

def 光轴(block, data_dir):

    # 光轴_所有点的角度偏差结果(mAA)
    df = pd.read_csv(data_dir/'c_mAA_r.csv')
    fig_args={
        'title':'光轴_所有点的角度偏差结果(mAA)',
    }
    fig = scatter_chart(df,**fig_args)
    block[0].plotly_chart(fig)
    # 光轴_所有点的位移偏差结果(normal)
    df = pd.read_csv(data_dir/'c_normal_t.csv')
    fig_args={
        'title':'光轴_所有点的位移偏差结果(normal)',
    }
    fig = scatter_chart(df,**fig_args)
    block[1].plotly_chart(fig)
    # 光轴_最后一个点的位移偏差结果
    # df = pd.read_csv(data_dir/'c_dif_per.csv')[['9']]
    df = pd.read_csv(data_dir/'c_normal_t.csv')[['9']]
    fig_args={
        'title':'光轴_最后一个点的位移偏差结果',
    }
    fig = scatter_chart(df,**fig_args)
    block[2].plotly_chart(fig)

def 焦距(block, data_dir):

    # 焦距_所有点的角度偏差结果(mAA)
    df = pd.read_csv(data_dir/'f_mAA_r.csv')
    fig_args={
        'title':'焦距_所有点的角度偏差结果(mAA)',
    }
    fig = scatter_chart(df,**fig_args)
    block[0].plotly_chart(fig)
    # 焦距_所有点的位移偏差结果(normal)
    df = pd.read_csv(data_dir/'f_normal_t.csv')
    fig_args={
        'title':'焦距_所有点的位移偏差结果(normal)',
    }
    fig = scatter_chart(df,**fig_args)    
    block[1].plotly_chart(fig)
    # 焦距_最后一个点的位移偏差结果
    # df = pd.read_csv(data_dir/'f_normal_t.csv')[['9']]
    df = pd.read_csv(data_dir/'f_normal_t.csv')[['9']]
    fig_args={
        'title':'焦距_最后一个点的位移偏差结果',
    }
    fig = scatter_chart(df,**fig_args)
    block[2].plotly_chart(fig)

def 特征点(block, data_dir):

    百分比 = int(st.slider('特征点偏移百分比', 5, 100, 5, 5))
    百分比 = '{:.1f}'.format(百分比*1.)
    
    # 所有点的角度偏差结果(mAA)
    df = list()
    for i in range(0,10):
        df.append(pd.read_excel(data_dir/'m_mAA_r.xlsx', sheet_name=i)[[百分比]].rename(columns={百分比: i}))
    fig_args={
        'title':'所有点的角度偏差结果(mAA)',
    }
    fig = scatter_chart(pd.concat(df, axis=1),**fig_args)
    block[0].plotly_chart(fig)
    # 所有点的位移偏差结果(normal)
    df = list()
    for i in range(0,10):
        df.append(pd.read_excel(data_dir/'m_normal_t.xlsx', sheet_name=i)[[百分比]].rename(columns={百分比: i}))
    print(pd.concat(df, axis=1))
    fig_args={
        'title':'所有点的位移偏差结果(normal)',
    }
    fig = scatter_chart(pd.concat(df, axis=1),**fig_args)
    block[1].plotly_chart(fig)
    # 最后一个点的位移偏差结果
    df = pd.read_csv(data_dir/'m_dif_per.csv')[[百分比]].rename(columns={百分比: i})
    fig_args={
        'title':'最后一个点的位移偏差结果',
    }
    fig = scatter_chart(df,**fig_args)
    block[2].plotly_chart(fig)

def display_plot(st, data_dir):
    data_dir = Path(data_dir)
    col1 = st.columns(3)
    光轴(col1, data_dir)

    col2 = st.columns(3)
    焦距(col2, data_dir)

    col3 = st.columns(3)
    特征点(col3, data_dir)

    
if __name__ == '__main__':
    display_plot(st, '/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_img/eg')