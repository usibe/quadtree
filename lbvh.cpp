#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <functional>
#include <algorithm>
#include <vector>
#include <stdlib.h>

using namespace std;

#define NMAX 1000000 // ノードの最大数
#define IMG_SIZE 512 // 画像サイズ
#define MAX_LEVEL 16 // ツリーの最大深さ (32bit Mortonキーでは各座標16bitまで)
#define MAX_PARTICLES (1 << 10) // 最大粒子数

typedef struct TreeNode {
    int id = -1;
    double cm[2] = {0.0, 0.0}; // center of mass
    double pos[2] = {0.0, 0.0}; // ノードの左下の座標
    double mass = 0.0;
    double size = 0.0; // ノードのサイズ（幅と高さは同じと仮定）
    int children[2] = {-1, -1}; // 子ノードのID（-1は子なしを表す）
} TreeNode;


typedef struct Particle{
    double pos[2];
    double mass;
} Particle, *pParticle;




uint32_t get_key(uint32_t x, uint32_t y){
    x = (x|(x<<8)) & 0x00ff00ff;
    x = (x|(x<<4)) & 0x0f0f0f0f;
    x = (x|(x<<2)) & 0x33333333;
    x = (x|(x<<1)) & 0x55555555;
    y = (y|(y<<8)) & 0x00ff00ff;
    y = (y|(y<<4)) & 0x0f0f0f0f;
    y = (y|(y<<2)) & 0x33333333;
    y = (y|(y<<1)) & 0x55555555;
    return (x | (y << 1));
}

int compactBits(int x){
    x &= 0x55555555;
    x = (x ^ (x >> 1)) & 0x33333333;
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}



void p2key(vector<Particle> &particle, vector<pair<uint32_t,int> > &key, int n){
    //各粒子についてmortonkeyを求める
    for(int i=0; i<n; i++){
        uint32_t x, y;
        uint32_t scale = 1u << MAX_LEVEL;
        x = (uint32_t)(scale * particle[i].pos[0]);
        y = (uint32_t)(scale * particle[i].pos[1]);
        key[i].first = get_key(x, y);
        key[i].second = i;
    }
    sort(key.begin(), key.end());
}

inline int lcp(const vector<pair<uint32_t,int>> &key, int i, int j){
    if(j < 0 || j >= (int)key.size()) return -1;
    uint32_t x = key[i].first ^ key[j].first;
    if(x == 0) return MAX_LEVEL * 2; // 全ビットが同じ場合は最大の共通接頭辞長を返す
    return __builtin_clz(x) - (32 - MAX_LEVEL * 2);  // 右詰め補正
}

int findSplit(const vector<pair<uint32_t,int>> &key, int left, int right){
    int split = left;
    int minLCP = INT32_MAX;

    // right は exclusive end。区間は [left, right)
    for(int i = left; i < right - 1; i++){
        int v = lcp(key, i, i+1);
        if(v < minLCP) {
            minLCP = v;
            split = i;
        }
    }
    return split + 1; // exclusive end を返す
}

// 二分探索で分割点を見つける。findSplitの効率化版
int findSpritBinary(const vector<pair<uint32_t,int>> &key, int left, int right){
    // leftとright-1のLCPを比較して、LCPが小さい方に分割点があることを利用して二分探索で分割点を見つける
    int lp = left, rp = right - 1;
    while(lp+1 < rp){
        int mid = (lp + rp) / 2;
        if(lcp(key, left, mid) < lcp(key, mid, right - 1)) {
            rp = mid;
        } else {
            lp = mid;
        }
    }
    return lp + 1; // exclusive end を返す
}


void calcNodePosSize(TreeNode *node, const vector<pair<uint32_t,int>> &key, int left, int right){
    // ノードのセルサイズの計算
    int common_prefix_len = lcp(key, left, right-1);
    int level = common_prefix_len / 2; // 1レベルごとに2ビットずつ分割されるため、レベルは共通接頭辞の長さの半分
    double cell_size = 1.0 / (1 << level); // レベルに応じたセルサイズの計算
    
    // ノードの位置の計算
    uint32_t morton_key = key[left].first;
    // buildLBVH内のpos計算を修正
    int shift = (MAX_LEVEL - level) * 2;          // 下位の「ツリー内位置」ビットを捨てる
    uint32_t cell_bits = (int)(morton_key >> shift);    // 上位level*2ビットが残る
    if(shift >= 32) cell_bits = 0; // シフトが32以上の場合は全ビットが捨てられるため、cell_bitsは0になる

    // prefix を x と y に分割して、ノードの左下の座標を計算
    int x = compactBits(cell_bits);
    int y = compactBits(cell_bits >> 1);

    node->size = cell_size;
    node->pos[0] = x * cell_size;
    node->pos[1] = y * cell_size;
}

void buildLBVH(vector<TreeNode> &nodes,
                vector<Particle> &particle,
                vector<pair<uint32_t,int>> &key,
                vector<vector<double>> &csum,
                int *nid, int left, int right)
{
    int id = *nid;
    nodes[id].id = id;
    // right は exclusive end。区間は [left, right)
    nodes[id].mass = csum[right][2] - csum[left][2];
    nodes[id].cm[0] = (csum[right][0] - csum[left][0]) / nodes[id].mass;
    nodes[id].cm[1] = (csum[right][1] - csum[left][1]) / nodes[id].mass;

    calcNodePosSize(&nodes[id], key, left, right);
    
    // 葉ノード
    if(right - left == 1) return;

    // デバッグ用に分割点とキーの範囲を表示
    printf("\n[Split] range=(%d,%d)\n", left, right);

    // 内部ノードの分割
    // int split = findSplit(key, left, right);
    int split = findSpritBinary(key, left, right);

    (*nid)++;
    nodes[id].children[0] = *nid;
    buildLBVH(nodes, particle, key, csum, nid, left, split);

    (*nid)++;
    nodes[id].children[1] = *nid;
    buildLBVH(nodes, particle, key, csum, nid, split, right);
}



void init_condition(vector<Particle> &particle, int n){
    for(int i=0; i<n; i++){
        particle[i].pos[0] = drand48();
        particle[i].pos[1] = drand48();
        particle[i].mass = 1.0/n;
    }
}

void init_csum(vector<Particle> &particle, vector<pair<uint32_t,int> > &key, vector<vector<double>> &csum, int n){
    for(int i=0; i<n; i++){
        int idx = key[i].second;
        csum[i+1][0] = csum[i][0] + particle[idx].pos[0] * particle[idx].mass;
        csum[i+1][1] = csum[i][1] + particle[idx].pos[1] * particle[idx].mass;
        csum[i+1][2] = csum[i][2] + particle[idx].mass;
    }
}




// PPM画像出力用関数群
void drawParticle(unsigned char img[IMG_SIZE][IMG_SIZE][3], double x, double y) {
    int px = (int)(x * IMG_SIZE);
    int py = (int)(y * IMG_SIZE);
    int radius = 4; // パーティクルの半径（ピクセル）
    for(int dy = -radius; dy <= radius; dy++) {
        for(int dx = -radius; dx <= radius; dx++) {
            int nx = px + dx;
            int ny = py + dy;
            if(nx >= 0 && nx < IMG_SIZE && ny >= 0 && ny < IMG_SIZE) {
                // 円形にする
                if(dx*dx + dy*dy <= radius*radius) {
                    img[ny][nx][0] = 255;
                    img[ny][nx][1] = 0;
                    img[ny][nx][2] = 0;
                }
            }
        }
    }
}

void drawBoundary(unsigned char img[IMG_SIZE][IMG_SIZE][3], double x, double y, double size) {
    int s = max(1, (int)(size * IMG_SIZE));
    int x0 = min(IMG_SIZE-1, (int)(x * IMG_SIZE));
    int y0 = min(IMG_SIZE-1, (int)(y * IMG_SIZE));
    int x1 = min(IMG_SIZE, x0 + s);
    int y1 = min(IMG_SIZE, y0 + s);
    // 上下
    for(int i = x0; i < x1; i++) {
        if(i >= 0 && i < IMG_SIZE && y0 >= 0 && y0 < IMG_SIZE) {
            img[y0][i][0] = 0;
            img[y0][i][1] = 255;
            img[y0][i][2] = 0;
        }
        if(i >= 0 && i < IMG_SIZE && y0 + s - 1 >= 0 && y0 + s - 1 < IMG_SIZE) {
            img[y0 + s - 1][i][0] = 0;
            img[y0 + s - 1][i][1] = 255;
            img[y0 + s - 1][i][2] = 0;
        }
    }
    // 左右
    for(int j = y0; j < y1; j++) {
        if(x0 >= 0 && x0 < IMG_SIZE && j >= 0 && j < IMG_SIZE) {
            img[j][x0][0] = 0;
            img[j][x0][1] = 255;
            img[j][x0][2] = 0;
        }
        if(x0 + s - 1 >= 0 && x0 + s - 1  < IMG_SIZE && j >= 0 && j < IMG_SIZE) {
            img[j][x0 + s - 1][0] = 0;
            img[j][x0 + s - 1][1] = 255;
            img[j][x0 + s - 1][2] = 0;
        }
    }
}

void drawTreeHierarchy(unsigned char img[IMG_SIZE][IMG_SIZE][3], vector<TreeNode> &nodes, int node_id){
    TreeNode *node = &nodes[node_id];

    drawBoundary(img, node->pos[0], node->pos[1], node->size);

    for(int i = 0; i < 2; i++){
        if(node->children[i] != -1){
            drawTreeHierarchy(img, nodes, node->children[i]);
        }
    }
}

void outputBitmap(const char *filename, vector<TreeNode> &node, int root_id, vector<Particle> &particle, int n) {
    unsigned char img[IMG_SIZE][IMG_SIZE][3];
    // 白で初期化
    for(int y = 0; y < IMG_SIZE; y++) {
        for(int x = 0; x < IMG_SIZE; x++) {
            img[y][x][0] = 255;
            img[y][x][1] = 255;
            img[y][x][2] = 255;
        }
    }
    // 境界線描画
    drawTreeHierarchy(img, node, root_id);
    // パーティクル描画
    for(int i = 0; i < n; i++) {
        drawParticle(img, particle[i].pos[0], particle[i].pos[1]);
    }
    // PPM出力
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "P3\n%d %d\n255\n", IMG_SIZE, IMG_SIZE);
    for(int y = 0; y < IMG_SIZE; y++) {
        for(int x = 0; x < IMG_SIZE; x++) {
            fprintf(fp, "%d %d %d ", img[y][x][0], img[y][x][1], img[y][x][2]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("Bitmap output to %s\n", filename);
}



int main(){
    int n;
    fprintf(stdout,"input number of particles > ");
    scanf("%d", &n);
    if(n <= 0){
        fprintf(stderr, "Error: the number of particles must be positive\n");
        exit(1);
    }
    if(n > MAX_PARTICLES){
        fprintf(stderr, "Error: the number of particles is too large\n");
        exit(1);
    }

    // left/right は exclusive end を使う
    // levelの数の階層までしか掘り下げられない.
    // levelはツリーの深さを決めるパラメータ。大きいほど細かく分割されるが、計算量も使用メモリも増える。
    // levelに応じて、最下層は4^(level-1)分割されることになる。
    int nid = 0;
    int left = 0;
    int right = n; // exclusive end

    vector<Particle> particle(n);
    vector<pair<uint32_t,int> > key(n);
    vector<TreeNode> nodes(NMAX);
    vector<vector<double>> csum(n+1, vector<double>(3, 0.0)); // 累積和配列(0: x*mass, 1: y*mass, 2: mass)
    
    init_condition(particle, n);
    p2key(particle, key, n);
    init_csum(particle, key, csum, n); 

    // ツリー構築
    buildLBVH(nodes, particle, key, csum, &nid, left, right);

    // 画像出力
    outputBitmap("output.ppm", nodes, 0, particle, n);
    return 0;
}
