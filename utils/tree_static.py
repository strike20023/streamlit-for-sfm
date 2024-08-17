# %%
import numpy as np
import hashlib
from graphviz import Digraph
import pandas as pd
from pathlib import Path
from itertools import product
# %%
def get_branch_data(branch_csv):
    col_df = pd.read_csv(branch_csv)
    branch_data = []
    for i in [t for t in [i[:-2] for i in col_df.columns if '_l' in i]]:
        branch_data.append([(a, b) for a, b in zip(
            col_df[i+'_l'].tolist(), col_df[i+'_r'].tolist())])
    branch_data = list(product(*branch_data))
    branch_data = [[t for k in i for t in k] for i in branch_data]
    branch_data = pd.DataFrame(branch_data, columns=col_df.columns)
    branch_data['接受'] = False
    branch_data['数量'] = ''
    branch_data['_数量'] = ''
    return branch_data
# %%
def get_product_data(branch_data, data_file):
    p = Path(data_file)
    df = pd.read_csv(p).drop(columns=['reward'])
    df['done'] = df['done'].apply(lambda x: 1 if x else 0)
    df
    # %%
    product_data = branch_data.copy()
    feature_cols = df.drop(columns=['done']).columns
    for idx, row in product_data.iterrows():
        count = df[
            (df[feature_cols[0]] >= row[feature_cols[0]+'_l']) &
            (df[feature_cols[0]] < row[feature_cols[0]+'_r']) &
            (df[feature_cols[1]] >= row[feature_cols[1]+'_l']) &
            (df[feature_cols[1]] < row[feature_cols[1]+'_r']) &
            (df[feature_cols[2]] >= row[feature_cols[2]+'_l']) &
            (df[feature_cols[2]] < row[feature_cols[2]+'_r']) &
            (df[feature_cols[3]] >= row[feature_cols[3]+'_l']) &
            (df[feature_cols[3]] < row[feature_cols[3]+'_r'])
        ]['done']
        refuse = count.sum()
        accept = count.count()-refuse
        product_data.loc[idx, '_数量'] = count.count()
        product_data.loc[idx, '数量'] = str([accept, refuse])
        if accept > refuse:
            product_data.loc[idx, '接受'] = True
    product_data = product_data[product_data['_数量'] > 0]
    # product_data = product_data[product_data['接受'] == True]
    return product_data
# %%


def tree_to_code(结果表):
    展示表 = list()
    保留特征 = [
        i.replace("_l", "")
        for i in 结果表.columns if i[-2:] == "_l"
    ]
    for _, row in 结果表.iterrows():
        字典 = dict()
        for 特征 in 保留特征:
            if ~np.isinf(row[特征 + "_l"]) and ~np.isinf(row[特征 + "_r"]):
                # 字典[特征] = f"{row[特征+'_l']} <= {特征} < {row[特征+'_r']}"
                字典[特征] = f"{特征}\n{row[特征+'_l']}~{row[特征+'_r']}"
            elif ~np.isinf(row[特征 + "_l"]):
                # 字典[特征] = f"{row[特征+'_l']} <= {特征}"
                字典[特征] = f"{特征}\n>{row[特征+'_l']}"
            elif ~np.isinf(row[特征 + "_r"]):
                # 字典[特征] = f"{特征} < {row[特征+'_r']}"
                字典[特征] = f"{特征}\n<{row[特征+'_r']}"
            else:
                字典[特征] = ""
        字典["接受"] = "接受" if row["接受"] else "拒绝"
        字典["数量"] = str(row["数量"])
        展示表.append(字典)
    展示表 = pd.DataFrame(展示表)
    return 展示表




def 创建树图(df, output_format):
    dot = Digraph(format=output_format)
    root = "分类树"
    dot.node(root, root, shape="box")
    level_nodes = {
        i: list() for i in range(len(df.columns))
    }
    edges = set()

    def add_edges(df, parent, level=0):

        if level >= len(df.columns):
            return
        if level == len(df.columns) - 1:
            for _, group in df.groupby(df.columns[level]):
                node = group[df.columns[level]].iloc[0]
                combined_str = df.to_string() + parent + str(level)
                nodeid = hashlib.sha256(combined_str.encode()).hexdigest()
                dot.node(nodeid, node, shape="box")
                level_nodes[level].append(nodeid)
                dot.edge(parent, nodeid)
                edges.add((parent, nodeid))
            return

        current_column = df.columns[level]
        next_column = df.columns[level + 1] if level + \
            1 < len(df.columns) else None

        for _, group in df.groupby(current_column):
            node = group[current_column].iloc[0]
            nodeid = parent + node
            if node == "":
                add_edges(group, parent, level + 1)
            else:
                if nodeid not in dot.body:
                    dot.node(nodeid, node, shape="box")
                    level_nodes[level].append(nodeid)
                if (parent, nodeid) not in edges:
                    dot.edge(parent, nodeid)
                    edges.add((parent, nodeid))
                if next_column:
                    add_edges(group, nodeid, level + 1)

    add_edges(df, root)

    for k, v in level_nodes.items():
        with dot.subgraph() as s:
            s.attr(rank='same')
            for i in v:
                s.node(i)
    return dot

def get_tree(branch_csv, data_file, output_format='pdf'):
    branch_data = get_branch_data(branch_csv)
    product_data = get_product_data(branch_data, data_file)
    展示表 = tree_to_code(product_data)
    dot = 创建树图(展示表, output_format)
    dot.render(Path(data_file).parent/'tree')
    return
if __name__ == '__main__':
    get_tree('分支.csv', '../result_tree/eg/results_1-5.csv', 'svg')