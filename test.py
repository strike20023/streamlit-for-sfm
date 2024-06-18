import streamlit as st
from streamlit import session_state as ss

st.title('hello')
ss.slider_value = 5
ss.btn = st.button('push')
if ss.btn:
    st.info('good!')
    slider_value = st.slider('slide', 0, 10, ss.slider_value)
    st.write(slider_value)

import streamlit as st

# 保存状态
# if st.button('Button 1'):
#     st.session_state.button1_clicked = True

# if st.button('Button 2'):
#     st.session_state.button2_clicked = True

# # 检索状态
# if 'button1_clicked' not in st.session_state:
#     st.session_state.button1_clicked = False
# if 'button2_clicked' not in st.session_state:
#     st.session_state.button2_clicked = False

# # 使用状态信息进行逻辑判断
# if st.session_state.button1_clicked:
#     st.write('Button 1 clicked')
# if st.session_state.button2_clicked:
#     st.write('Button 2 clicked')
def check_button_clicked(button_name):
    if st.button(button_name):
        st.session_state[button_name] = True
    if button_name not in st.session_state:
        st.session_state[button_name] = False
    return st.session_state[button_name]
sss = st.checkbox('sss')
if sss:
    st.write(sss)
    if check_button_clicked('Button 1'):
        st.write('Button 1 clicked')
        import time
        st.write(time.time())
        s = st.slider('slider', 0, 10, 5)
        st.write(s)
if check_button_clicked('Button 2'):
    st.write('Button 2 clicked')

with st.form('my_form'):
    button1_clicked = st.form_submit_button('Button 1')
    button2_clicked = st.form_submit_button('Button 2')

if button1_clicked:
    st.write('Button 1 clicked')

if button2_clicked:
    st.write('Button 2 clicked')
