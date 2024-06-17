# %%
import os
import numpy as np
from treelib import Tree
from graphviz import Digraph
import pandas as pd

class SetManager:
    def __init__(self, interval):
        self.interval = interval
        self.sets = []

    def __add__(self, num):
        merged = False
        for s in self.sets:
            if any(abs(num - x) <= self.interval for x in s):
                s.add(num)
                merged = True
                break
        if not merged:
            self.sets.append({num})
        self.merge_sets()
        return self

    def merge_sets(self):
        merged = True
        while merged:
            merged = False
            for i in range(len(self.sets)):
                for j in range(i + 1, len(self.sets)):
                    if any(abs(x - y) <= self.interval for x in self.sets[i] for y in self.sets[j]):
                        self.sets[i].update(self.sets[j])
                        self.sets.pop(j)
                        merged = True
                        break
                if merged:
                    break

    def __repr__(self):
        return f"{self.sets}"

# 使用graphviz可视化
def create_graph(tree):
    # colors = ["#FFD700", "#FF6347", "#4682B4", "#32CD32", "#8A2BE2", "#FFA500"]
    # colors = ["#F8BBD0", "#E1BEE7", "#BBDEFB", "#C8E6C9", "#FFF9C4", "#FFE0B2"]
    colors = ["#FAD9E2", "#D7BBE8", "#A4CDEB", "#CCF2D2", "#FFE5B4"]
    level_name = ["root", "variance", "percentage", "fxy", "cxy"]
    dot = Digraph(
        comment="My Tree",
        graph_attr={
            "rankdir": "LR",
        },
    )
    for node in tree.all_nodes_itr():
        level = tree.level(node.identifier)  # 获取节点的层级
        # if node.is_leaf():
        #     color = colors[-1]
        # else:
        color = colors[level % len(colors)]  # 根据层级选择颜色
        if node.is_root():
            label = node.tag
        else:
            label = f'<<FONT POINT-SIZE="7"> {level_name[level]}</FONT><BR/><FONT POINT-SIZE="10">{node.tag}</FONT>>'
        dot.node(
            node.identifier, label=label, shape="box", style="filled", fillcolor=color
        )
        if not node.is_root():
            dot.edge(tree.parent(node.identifier).identifier, node.identifier)

    dot.render("result_trees/tree-test.dot", format="svg")


def split_param(df, idx_p):
    param = params[idx_p]
    if df["done"].sum() == 0:
        return df
    p_list = df[~df["done"]][param].sort_values().tolist()
    if idx_p != len(params) - 1:
        split_p = np.round(np.linspace(0, len(p_list) - 1, param_num[idx_p] + 1), 0)
        p_list = [p_list[int(i)] for i in split_p]
        for i in range(len(p_list) - 1):
            if i == len(p_list) - 2:
                res = df[(df[param] >= p_list[i]) & (df[param] <= p_list[i + 1])]
            else:
                res = df[(df[param] >= p_list[i]) & (df[param] < p_list[i + 1])]
            if len(res) > 0:
                yield res
    else:
        sm = SetManager(interval=minimum_interval[idx_p])
        for i in p_list:
            sm += i
        max_sets = sorted(sm.sets, key=len, reverse=True)
        # if len(max_sets)>1:
        #     print(max_sets)
        #     for i in max_sets: 
        #         print(df[df[param].isin(i)])
        for i in max_sets: 
            yield df[df[param].isin(i)]
# %%

# %%

def get_label(df, idx_p):
    min_max = df[~df["done"]][params[idx_p]].describe()[["min", "max"]].tolist()
    if min_max[0] != min_max[1]:
        return "{}&rarr;{}".format(*min_max)
    else:
        return min_max[0]


def create_nodes(df):
    tree = Tree()
    tree.create_node("root", "root")
    node_name_format = "{}_{}".format

    # try:
    with pd.option_context("display.max_rows", None):
        for idx_i, i in enumerate(split_param(df, 0)):
            label = get_label(i, 0)
            level1 = tree.create_node(
                label, node_name_format("root", idx_i), parent="root"
            ).identifier
            sub_df = [sub for sub in split_param(i, 1)]
            for idx_j, j in enumerate(sub_df):
                label = get_label(j, 1)
                level2 = tree.create_node(
                    label, node_name_format(level1, idx_j), parent=level1
                ).identifier
                sub_df = [sub for sub in split_param(j, 2)]
                for idx_k, k in enumerate(sub_df):
                    label = get_label(k, 2)
                    sub_df = [sub for sub in split_param(k, 3)]
                    if sub_df:
                        level3 = tree.create_node(
                            label, node_name_format(level2, idx_k), parent=level2
                        ).identifier
                        for idx_l, l in enumerate(sub_df):
                            label = get_label(l, 3)
                            tree.create_node(
                                label, node_name_format(level3, idx_l), parent=level3
                            )

    return tree


# 创建一个树对象
def create_tree(csv_file, branch_num):
    global params, param_num, minimum_interval
    param_num = branch_num
    params = ["variance", "percentage", "fxy", "cxy"]
    os.getcwd()
    df = pd.read_csv(csv_file)
    get_minimum_interval = lambda x:round(min([abs(x[i]-x[i-1]) for i,_ in enumerate(range(1,len(x)))]), 4)
    minimum_interval = [get_minimum_interval(df[param].sort_values().drop_duplicates().tolist()) for param in params]
    tree = create_nodes(df)
    create_graph(tree)


if __name__ == "__main__":
    create_tree("/mnt/ST8000/zhenhanbai/RL/Results/results_2.csv", [4,4,4,4])
