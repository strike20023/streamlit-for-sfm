# %%
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - {%(message)s}')
logging.info("123")
print('adwverve')
# %%
import argparse
args = argparse.ArgumentParser()
# %%
args.add_argument('--input_img', default=None, type=str, nargs='+', help="图片组编码/地址")