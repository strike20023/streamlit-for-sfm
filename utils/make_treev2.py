# %%
import hashlib
from graphviz import Digraph
from sklearn.tree import _tree
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
import graphviz
import numpy as np
def 绘制决策树(csv_file,dot_path):
    print(dot_path.replace('.dot','.2.dot'))
    数据集 = pd.read_csv(csv_file).drop(columns=["reward"])

    特征 = 数据集.drop(columns=["done"])
    标签 = ~数据集["done"]

    分类器 = DecisionTreeClassifier(max_depth=6)
    分类器 = 分类器.fit(特征, 标签)

    dot数据 = tree.export_graphviz(
        分类器,
        out_file=None,
        feature_names=特征.columns,
        class_names=["否", "是"],
        filled=True,
        rounded=True,
        special_characters=True,
        leaves_parallel=True,
    )
    图形 = graphviz.Source(dot数据)
    图形.render(dot_path.replace('.dot','.2.dot'))
    print(dot_path.replace('.dot','.2.dot'))

    def tree_to_code(tree, feature_names):
        tree_ = tree.tree_
        feature_name = [
            feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
            for i in tree_.feature
        ]
        paths = []

        def recurse(node, depth, path):
            indent = "  " * depth
            if tree_.feature[node] != _tree.TREE_UNDEFINED:
                name = feature_name[node]
                threshold = tree_.threshold[node]
                path_left = [name, "<=", threshold]
                path_right = [name, ">", threshold]
                recurse(tree_.children_left[node], depth + 1, path + [path_left])
                recurse(tree_.children_right[node], depth + 1, path + [path_right])
            else:
                paths.append((path, tree_.value[node]))

        recurse(0, 0, [])
        return paths


    paths = tree_to_code(分类器, 特征.columns)
    结果表 = pd.DataFrame(
        columns=[*[t for i in 特征.columns for t in [i + "_l", i + "_r"]], "接受"]
    )
    for path, value in paths:
        row = dict()
        for node in path:
            node[-1] = round(node[-1], 3)
            # print(node, end=" ")
            suffix = "_r" if node[1] == "<=" else "_l"
            if node[0] + suffix in row:
                if node[1] == "<=":
                    row[node[0] + suffix] = min(row[node[0] + suffix], node[-1])
                else:
                    row[node[0] + suffix] = max(row[node[0] + suffix], node[-1])
            else:
                row[node[0] + suffix] = node[-1]
        row["接受"] = value[0][0] == 0
        结果表.loc[len(结果表)] = row
    结果表 = 结果表.dropna(axis=1, how="all")
    display(结果表)

    展示表 = list()
    保留特征 = set(
        [
            i.replace("_l", "").replace("_r", "")
            for i in 结果表.drop(columns=["接受"]).columns
        ]
    )
    for _, row in 结果表.iterrows():
        字典 = dict()
        for 特征 in 保留特征:
            if ~np.isnan(row[特征 + "_l"]) and ~np.isnan(row[特征 + "_r"]):
                # 字典[特征] = f"{row[特征+'_l']} <= {特征} < {row[特征+'_r']}"
                字典[特征] = f"{特征}\n{row[特征+'_l']}~{row[特征+'_r']}"
            elif ~np.isnan(row[特征 + "_l"]):
                # 字典[特征] = f"{row[特征+'_l']} <= {特征}"
                字典[特征] = f"{特征}\n>={row[特征+'_l']}"
            elif ~np.isnan(row[特征 + "_r"]):
                # 字典[特征] = f"{特征} < {row[特征+'_r']}"
                字典[特征] = f"{特征}\n<{row[特征+'_r']}"
            else:
                字典[特征] = ""
        字典["接受"] = "接受" if row["接受"] else "拒绝"
        展示表.append(字典)
    展示表 = pd.DataFrame(展示表)
    展示表 = 展示表


    def 创建树图(df):
        dot = Digraph(format="svg")
        root = "Root"
        dot.node(root, root, shape="box")
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
                    dot.edge(parent, nodeid)
                    edges.add((parent, nodeid))
                return

            current_column = df.columns[level]
            next_column = df.columns[level + 1] if level + 1 < len(df.columns) else None

            for _, group in df.groupby(current_column):
                node = group[current_column].iloc[0]
                nodeid = parent + node
                if node == "":
                    add_edges(group, parent, level + 1)
                else:
                    if nodeid not in dot.body:
                        dot.node(nodeid, node, shape="box")
                    if (parent, nodeid) not in edges:
                        dot.edge(parent, nodeid)
                        edges.add((parent, nodeid))
                    if next_column:
                        add_edges(group, nodeid, level + 1)

        add_edges(df, root)
        return dot


    dot = 创建树图(展示表)
    dot.render(dot_path)
if __name__ == '__main__':
    from pathlib import Path
    for p in Path('/mnt/ST8000/huangzhe/海纹项目/streamlit_app/result_trees/').glob('*.csv'):
        绘制决策树(p,p.__str__().replace('.csv','.dot'))
# %%
