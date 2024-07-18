# %%
import numpy as np
import hashlib
from graphviz import Digraph
import pandas as pd
from pathlib import Path


class DecisionTreeNode:
    def __init__(self, gini, num_samples, num_samples_per_class, predicted_class, feature_ranges=None):
        self.gini = gini
        self.num_samples = num_samples
        self.num_samples_per_class = num_samples_per_class
        self.predicted_class = predicted_class
        self.feature_index = None
        self.threshold = None
        self.left = None
        self.right = None
        self.feature_ranges = feature_ranges if feature_ranges else []

    def __repr__(self):
        return (f"DecisionTreeNode({self.predicted_class}, num_samples_per_class={self.num_samples_per_class}, gini={self.gini}, samples={self.num_samples}, "
                f"feature_index={self.feature_index}, threshold={self.threshold}), feature_ranges={self.feature_ranges}")

    def to_dot(self, graph=[], parent_id=None, edge_label=""):
        node_id = len(graph)
        # 生成当前节点的标签
        labels = [f"{cls}: {count}" for cls, count in zip(
            ["True", "False"], self.num_samples_per_class)]
        # feature_ranges
        labels.append(f"{self.feature_ranges}")
        label = f"n{node_id} [label=\"{'\\n'.join(labels)}\"]"
        graph.append(label)
        if parent_id is not None:
            # 为当前节点和父节点之间添加边
            graph.append(
                f"n{parent_id} -> n{node_id} [label=\"{edge_label}\"]")

        # 递归为子节点生成DOT描述
        if self.left:
            self.left.to_dot(graph, node_id, f"{
                             self.feature_index} <= {self.threshold}")
        if self.right:
            self.right.to_dot(graph, node_id, f"{
                              self.feature_index} > {self.threshold}")
        return graph


def gini_impurity(labels):
    """计算基尼不纯度"""
    classes = set(labels)
    impurity = 1
    for cls in classes:
        prob_cls = labels.count(cls) / len(labels)
        impurity -= prob_cls ** 2
    return impurity


def split(dataset, feature_index, threshold):
    """根据特征索引和阈值分割数据集"""
    left = [sample for sample in dataset if sample[0]
            [feature_index] <= threshold]
    right = [sample for sample in dataset if sample[0]
             [feature_index] > threshold]
    return left, right


def best_split(dataset, feature_indexs):
    """寻找最佳分割"""
    best_index, best_threshold, best_gini, best_splits = None, None, float(
        'inf'), None
    for feature_index in feature_indexs:  # 遍历每个特征
        thresholds = set([sample[0][feature_index] for sample in dataset])
        for threshold in thresholds:
            left, right = split(dataset, feature_index, threshold)
            if not left or not right:
                continue
            gini = (len(left) / len(dataset) * gini_impurity([x[1] for x in left]) +
                    len(right) / len(dataset) * gini_impurity([x[1] for x in right]))
            if gini < best_gini:
                best_gini, best_index, best_threshold, best_splits = gini, feature_index, threshold, (
                    left, right)
    return best_index, best_threshold, best_splits


def build_tree(dataset, num_selected_features, depth=0, max_depth=1000, feature_ranges=None):
    """递归构建决策树"""
    if feature_ranges is None:
        feature_ranges = [(float('inf'), float('-inf'))] * \
            len(num_selected_features)

    labels = [sample[1] for sample in dataset]
    num_samples_per_class = [labels.count(i) for i in [True, False]]
    predicted_class = max(set(labels), key=labels.count)
    node = DecisionTreeNode(
        gini=gini_impurity(labels),
        num_samples=len(dataset),
        num_samples_per_class=num_samples_per_class,
        predicted_class=predicted_class,
        feature_ranges=feature_ranges
    )
    # 递归终止条件

    ratio = 0.1

    rate = num_samples_per_class[0] / \
        (num_samples_per_class[0]+num_samples_per_class[1])
    if rate < ratio or rate > (1-ratio):
        return node
    if len(labels) < 100:
        return node

    if depth >= max_depth:
        print('max_depth', depth)
        return node
    if not all(num_samples_per_class):
        print('num_samples_per_class', num_samples_per_class)
        return node
    feature_index = num_selected_features[0:1]
    index, threshold, splits = best_split(dataset, feature_index)
    if splits:
        pass
    elif all(num_samples_per_class):
        # num_selected_features = num_selected_features[1:]
        # feature_index = num_selected_features[0:1]
        # index, threshold, splits = best_split(dataset, feature_index)
        print(num_selected_features)
        if num_selected_features:
            return build_tree(dataset, num_selected_features[1:], depth + 1, max_depth, feature_ranges)
        else:
            return node
    else:
        return node
    node.feature_index = index
    node.threshold = threshold
    left, right = splits
    print('+', num_samples_per_class)
    l_feature_ranges = feature_ranges.copy()
    r_feature_ranges = feature_ranges.copy()
    l_feature_ranges[index] = (feature_ranges[index][0], threshold)
    r_feature_ranges[index] = (threshold, feature_ranges[index][1])
    node.left = build_tree(left, num_selected_features,
                           depth + 1, max_depth, l_feature_ranges)
    node.right = build_tree(right, num_selected_features,
                            depth + 1, max_depth, r_feature_ranges)
    return node


def predict(sample, tree):
    """使用决策树预测单个样本"""
    while tree.left:
        if sample[tree.feature_index] <= tree.threshold:
            tree = tree.left
        else:
            tree = tree.right
    return tree.predicted_class


# dataset = [(np.random.rand(3), np.random.choice([True, False])) for _ in range(6)]
# from pprint import pprint
# tree = build_tree(dataset, num_selected_features=[0 for i in range(len(dataset[0][0]))])


def print_dot(tree):
    # 初始化DOT描述
    dot_graph = ['digraph Tree {', 'node [shape=box];']
    # 添加节点和边
    dot_graph.extend(tree.to_dot())
    dot_graph.append('}')
    # 打印或返回DOT格式的字符串
    return "\n".join(dot_graph)


def print_leaf_nodes(node):
    leaf_info = []
    if node is not None:
        if node.left is None and node.right is None:  # 叶子节点
            # print(node)
            leaf_info.append(
                (node.feature_ranges, node.predicted_class, node.num_samples_per_class))
        else:
            leaf_info.extend(print_leaf_nodes(node.left))
            leaf_info.extend(print_leaf_nodes(node.right))
    return leaf_info
# dot_representation = print_dot(tree)


# %%
# p = Path('../result_tree/exp/2024-06-07-19:30/results_6-5.csv')
p = Path('/Users/huangzhe/Documents/GitHub/streamlit-for-sfm/result_tree/eg/results_1-5.csv')
# p = Path('/Users/huangzhe/Documents/GitHub/streamlit-for-sfm/result_tree/eg/results.csv')

数据集 = pd.read_csv(p).drop(columns=["reward"])

特征 = 数据集.drop(columns=["done"])[['fxy', 'cxy', 'percentage', 'variance']]
# [['fxy', 'percentage', 'variance', 'cxy']]
标签 = ~数据集["done"]
dataset = [(特征.to_numpy()[i].tolist(), 标签.to_numpy()[i])
           for i in range(len(特征))]


# pprint(dataset)
tree = build_tree(dataset, num_selected_features=[
                  i for i in range(len(dataset[0][0]))])
# tree = build_tree(dataset, num_selected_features=[
#                   i for i in range(len(dataset[0][0]))][::-1])
# tree = build_tree(dataset, num_selected_features=order)
dot_representation = print_dot(tree)
# print(dot_representation)
结果表 = print_leaf_nodes(tree)
print(len(结果表))
结果表 = [[*[k for t in i[0] for k in t], i[1], i[2]] for i in 结果表]
结果表 = pd.DataFrame(
    结果表, columns=[*[i+t for i in 特征.columns for t in ['_l', '_r']], "接受", "数量"])
# 结果表 = 结果表[结果表["接受"] == True]
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


展示表 = tree_to_code(结果表)


def 创建树图(df):
    dot = Digraph(format="pdf")
    root = "Root"
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


dot = 创建树图(展示表)
dot.render('./tree')
# %%
展示表
# %%
