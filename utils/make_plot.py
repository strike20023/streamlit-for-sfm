import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

st.markdown("""<style>#MainMenu {visibility: hidden;} footer {visibility: hidden;}</style>""", unsafe_allow_html=True)

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
    st.plotly_chart(fig)


def display_plot(st, data_dir):
    data_dir = Path(data_dir)

    # 光轴_所有点的角度偏差结果(mAA)
    df = pd.read_csv(data_dir/'c_mAA_r.csv')
    fig_args={
        'title':'光轴_所有点的角度偏差结果(mAA)',
    }
    scatter_chart(df,**fig_args)
    # 光轴_所有点的位移偏差结果(normal)
    df = pd.read_csv(data_dir/'c_normal_t.csv')
    fig_args={
        'title':'光轴_所有点的位移偏差结果(normal)',
    }
    scatter_chart(df,**fig_args)
    # 光轴_最后一个点的位移偏差结果
    # df = pd.read_csv(data_dir/'c_dif_per.csv')[['9']]
    df = pd.read_csv(data_dir/'c_normal_t.csv')[['9']]
    fig_args={
        'title':'光轴_最后一个点的位移偏差结果',
    }
    scatter_chart(df,**fig_args)

    # 焦距_所有点的角度偏差结果(mAA)
    df = pd.read_csv(data_dir/'f_mAA_r.csv')
    fig_args={
        'title':'焦距_所有点的角度偏差结果(mAA)',
    }
    scatter_chart(df,**fig_args)
    # 焦距_所有点的位移偏差结果(normal)
    df = pd.read_csv(data_dir/'f_normal_t.csv')
    fig_args={
        'title':'焦距_所有点的位移偏差结果(normal)',
    }
    scatter_chart(df,**fig_args)    
    # 焦距_最后一个点的位移偏差结果
    # df = pd.read_csv(data_dir/'f_normal_t.csv')[['9']]
    df = pd.read_csv(data_dir/'f_normal_t.csv')[['9']]
    fig_args={
        'title':'焦距_最后一个点的位移偏差结果',
    }
    scatter_chart(df,**fig_args)

    百分比 = int(st.slider('特征点偏移百分比', 5, 100, 5, 5))
    百分比 = '{:.1f}'.format(百分比*1.)
    print(百分比)

    # 所有点的角度偏差结果(mAA)
    df = list()
    for i in range(0,10):
        df.append(pd.read_excel(data_dir/'m_mAA_r.xlsx', sheet_name=i)[[百分比]].rename(columns={百分比: i}))
    fig_args={
        'title':'所有点的角度偏差结果(mAA)',
    }
    scatter_chart(pd.concat(df, axis=1),**fig_args)
    # 所有点的位移偏差结果(normal)
    df = list()
    for i in range(0,10):
        df.append(pd.read_excel(data_dir/'m_normal_t.xlsx', sheet_name=i)[[百分比]].rename(columns={百分比: i}))
    print(pd.concat(df, axis=1))
    fig_args={
        'title':'所有点的位移偏差结果(normal)',
    }
    scatter_chart(pd.concat(df, axis=1),**fig_args)
    # 最后一个点的位移偏差结果
    df = pd.read_csv(data_dir/'m_dif_per.csv')[[百分比]].rename(columns={百分比: i})
    fig_args={
        'title':'最后一个点的位移偏差结果',
    }
    scatter_chart(df,**fig_args)
if __name__ == '__main__':
    display_plot(st, '/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_img/eg')