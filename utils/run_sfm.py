# %% [markdown]
# ### 单因素分析实验

# %%
import os
import sys
import time
import shutil
import argparse
import pandas as pd
import logging as log

# %%
# python /mnt/ST8000/jialei/SFM_E/API/run_sfm.py --input_img [字符串] --variance [数字] --percentage [数字] --fxy [数字] --сху [数字] --output_excel_path [字符串]  --output_img_dir [字符串]
# python /mnt/ST8000/jialei/SFM_E/API/run_sfm.py  --output_excel_path /mnt/ST8000/jialei/SFM_E/API/result.xlsx --output_img_dir /mnt/ST8000/jialei/SFM_E/API/pic

# %%
def play_SfM(args):
    root_data_path = r"/mnt/ST8000/jialei/SFM_E/ex_sea_11/8-8"
    match_path = os.path.join(root_data_path, "match_points/output")
    focal_path = os.path.join(root_data_path, "focal_length/output")
    cxy_path = os.path.join(root_data_path, "cxcy/output")

    df_match = pd.read_csv(os.path.join(match_path , "dif_per.csv")
                         , header=0
                         , encoding="utf_8_sig")
    df_focal = pd.read_csv(os.path.join(focal_path , "dif_per.csv")
                         , header=0
                         , encoding="utf_8_sig")
    df_cxy = pd.read_csv(os.path.join(cxy_path , "dif_per.csv")
                         , header=0
                         , encoding="utf_8_sig")
    
    
    ew_sam = pd.ExcelWriter(args.output_excel_path)
    df_match.to_excel(ew_sam, sheet_name="match_points", index=False, header=True)
    df_focal.to_excel(ew_sam, sheet_name="focal_length", index=False, header=True)
    df_cxy.to_excel(ew_sam, sheet_name="cxcy", index=False, header=True)

    ew_sam.save()


# %%
def analysic_get_pic(args):
    # 分析数据
    pic_path = r"/mnt/ST8000/jialei/SFM_E/API/pic"
    target_pic_path = args.output_img_dir
    if not os.path.exists(target_pic_path):
        os.mkdir(target_pic_path)

    for file_name in os.listdir(pic_path):
        curr_path = os.path.join(pic_path, file_name)
        if file_name.endswith('.svg') or file_name.endswith('.png') or file_name.endswith('.jpg'):
            shutil.copy(curr_path, os.path.join(target_pic_path, file_name))


# %%
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_img', default=None, type=str, help="图片组编码/地址")
    parser.add_argument('--variance', default=None,type=float, nargs='+', help="方差")
    parser.add_argument('--percentage', default=None,type=float, nargs='+',help="匹配特征点误差分布百分比")
    parser.add_argument('--fxy', default=None,type=float, nargs='+',help="焦距")
    parser.add_argument('--cxy', default=None,type=float, nargs='+',help="光轴")
    parser.add_argument('--output_excel_path', default=None, type=str, help="输出excel地址")
    parser.add_argument('--output_img_dir', default=None, type=str, help="输出图片地址")
    args = parser.parse_args()
    # args.output_excel_path = r"/mnt/ST8000/jialei/SFM_E/API/result.xlsx"
    # args.output_img_dir = r"/mnt/ST8000/jialei/SFM_E/API/pic"
    
    # print(args)
    log.basicConfig(level=log.INFO, format='输出: {%(message)s}')
    
    for i in range(1,101):
        log.info(f"{i}%")


    play_SfM(args)

    analysic_get_pic(args)




# %%



