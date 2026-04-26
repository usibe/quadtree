# quadtree

## 概要
このプロジェクトは2次元粒子系のクワッドツリー構築アルゴリズムのサンプルです。粒子分布をツリー構造で管理し、画像として出力します。

### ファイル構成
- quadtree_basic.cpp: 通常の方法によるクワッドツリー構築
- quadtree_morton.cpp: Mortonキー（Zオーダー）を用いたクワッドツリー構築
- lbvh.cpp: Mortonキーを用いたMortonキーを利用したLBVH（Linear Bounding Volume Hierarchy）構築

## 使い方
1. コンパイル

```
g++ quadtree_basic.cpp -o quadtree_basic
```
または
```
g++ quadtree_morton.cpp -o quadtree_morton
```
または
```
g++ lbvh.cpp -o lbvh
```

2. 実行

```
./quadtree_basic
```
または
```
./quadtree_morton
```
または
```
./lbvh
```

粒子数を入力すると、output.ppm という画像ファイルが生成されます。

## 参考
- クワッドツリー
- Mortonキー（Zオーダー）
- LBVH（Linear Bounding Volume Hierarchy）

## ライセンス
MIT License
