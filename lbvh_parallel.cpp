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
#include <omp.h>

using namespace std;

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


pair<int, int> determineRange(const vector<pair<uint32_t,int>> &key, int idx, int n) {
    int lcp_left = lcp(key, idx, idx - 1);
    int lcp_right = lcp(key, idx, idx + 1);
    int direction = (lcp_left > lcp_right) ? -1 : 1; // 左右どちらに範囲が広がるか
    int lcp_min = (direction == 1) ? lcp_left : lcp_right; // 範囲の広がりを決定するLCPの最小値

    // 担当範囲の長さを二分探索で決める
    int lmax = 2;
    while(true){
        int next_idx = idx + direction * lmax;
        if(next_idx < 0 || next_idx >= n || lcp(key, idx, next_idx) <= lcp_min) break;
        lmax *= 2;
    }

    // 二分探索で正確な範囲を見つける
    int l = 0;
    for(int step = lmax/2; step > 0; step /= 2) {
        int next_idx = idx + direction * (l + step);
        if(next_idx >= 0 && next_idx < n && lcp(key, idx, next_idx) > lcp_min) {
            l += step;
        }
    }
    int left = min(idx, idx + direction * l);
    int right = max(idx, idx + direction * l) + 1; // exclusive end

    return {left, right}; // [left, right) の範囲を返す
}


void buildLBVH_parallel(vector<TreeNode> &nodes,
                vector<Particle> &particle,
                vector<pair<uint32_t,int>> &key, int n)
{
    #pragma omp parallel for
    for(int i=0; i<n; i++){ // 葉ノードの初期化
        int pid = key[i].second;
        nodes[n-1+i].id = n-1+i;
        nodes[n-1+i].mass = particle[pid].mass;
        nodes[n-1+i].cm[0] = particle[pid].pos[0];
        nodes[n-1+i].cm[1] = particle[pid].pos[1];
        nodes[n-1+i].pos[0] = particle[pid].pos[0];
        nodes[n-1+i].pos[1] = particle[pid].pos[1];
        nodes[n-1+i].size = 0.0; // 葉ノードのサイズは0とする
        nodes[n-1+i].children[0] = -1; // 葉ノードは子なし
        nodes[n-1+i].children[1] = -1;
    }

    #pragma omp parallel for
    for(int i = 0; i < n-1; i++) { // 枝ノード数は粒子数-1
        auto [left, right] = determineRange(key, i, n);
        int split = findSpritBinary(key, left, right);
        calcNodePosSize(&nodes[i], key, left, right);

        // 枝ノードはsplitを子ノードの分割点とする。splitはexclusive endなので、split-1が左の子の最後の葉ノードになる。
        // なぜ分岐点を左右の子ノードの境界にするのか？ -> これにより、分割された空間が重ならず、効率的なツリー構造が得られるため。
        if(split - left == 1){
            nodes[i].children[0] = (n-1)+left; // 左の子は葉ノード
        }else{
            nodes[i].children[0] = split-1; // 左の子は枝ノード
        }

        if(right - split == 1){
            nodes[i].children[1] = (n-1)+(right-1); // 右の子は葉ノード
        }else{
            nodes[i].children[1] = split; // 右の子は枝ノード
        }
    }

    // 枝ノードの位置とサイズの計算
    // 並列化は難しいため、ここはシングルスレッドで処理する
    for(int i = n-2; i >= 0; i--) {
        int left = nodes[i].children[0];
        int right = nodes[i].children[1];
        nodes[i].mass = nodes[left].mass + nodes[right].mass;
        nodes[i].cm[0] = (nodes[left].cm[0] * nodes[left].mass + nodes[right].cm[0] * nodes[right].mass) / nodes[i].mass;
        nodes[i].cm[1] = (nodes[left].cm[1] * nodes[left].mass + nodes[right].cm[1] * nodes[right].mass) / nodes[i].mass;
    }
}



void init_condition(vector<Particle> &particle, int n){
    for(int i=0; i<n; i++){
        particle[i].pos[0] = drand48();
        particle[i].pos[1] = drand48();
        particle[i].mass = 1.0/n;
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

    // OpenMPのテストコード
    // #pragma omp parallel for
    // for (int i = 0; i < 8; i++) {
    //     // 現在のスレッド番号を取得
    //     int tid = omp_get_thread_num();
    //     printf("Iteration %d is executed by thread %d\n", i, tid);
    // }

    vector<Particle> particle(n);
    vector<pair<uint32_t,int> > key(n);
    vector<TreeNode> nodes(2*n-1);    
    init_condition(particle, n);
    p2key(particle, key, n);

    // ツリー構築
    buildLBVH_parallel(nodes, particle, key, n);

    // 画像出力
    outputBitmap("output.ppm", nodes, 0, particle, n);
    return 0;
}
