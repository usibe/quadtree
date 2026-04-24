#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <functional>
#include <algorithm>
#include <vector>

using namespace std;

#define NMAX 1000000 // ノードの最大数
#define IMG_SIZE 512 // 画像サイズ
#define MAX_LEVEL 10 // ツリーの最大深さ
#define MAX_PARTICLES (1 << 10) // 最大粒子数

typedef struct TreeNode {
    int id = -1;
    double cm[2] = {0.0, 0.0}; // center of mass
    double pos[2] = {0.0, 0.0}; // ノードの左下の座標
    double mass = 0.0;
    double size = 0.0; // ノードのサイズ（幅と高さは同じと仮定）
    int children[4] = {-1, -1, -1, -1}; // 子ノードのID（-1は子なしを表す）
} TreeNode;


typedef struct Particle{
    double pos[2];
    double mass;
} Particle, *pParticle;




int get_key(int x, int y, int level){
    x = (x|(x<<8)) & 0x00ff00ff;
    x = (x|(x<<4)) & 0x0f0f0f0f;
    x = (x|(x<<2)) & 0x33333333;
    x = (x|(x<<1)) & 0x55555555;
    y = (y|(y<<8)) & 0x00ff00ff;
    y = (y|(y<<4)) & 0x0f0f0f0f;
    y = (y|(y<<2)) & 0x33333333;
    y = (y|(y<<1)) & 0x55555555;
    return (x|(y<<1));
}



void p2key(vector<Particle> &particle, vector<pair<int,int> > &key, int n, int level){
    //各粒子についてmortonkeyを求める
    for(int i=0; i<n; i++){
        int x, y, scale=1<<level; //levelの扱いがバラバラになってるので見直しが必要
        x = scale * particle[i].pos[0];
        y = scale * particle[i].pos[1];
        key[i].first = get_key(x, y, level);
        key[i].second = i;
    }
    sort(key.begin(), key.end());
}


void treeConstruction(vector<TreeNode> &node,
                     vector<Particle> &particle,
                     vector<pair<int,int> > &key, 
                     vector<vector<double>> &csum,
                     int *nid, int level, int left, int right){
    
    level -= 1;
    if(level < 0) return;
    int id = *nid;
    
    node[id].mass = csum[right][2] - csum[left][2];
    node[id].cm[0] = (csum[right][0] - csum[left][0]) / node[id].mass;
    node[id].cm[1] = (csum[right][1] - csum[left][1]) / node[id].mass;
    node[id].id = id; //いるかこれ？

    int start = left;
    int prev = (key[left].first >> (2*level)) & 0x3;
    for(int i=left+1; i<=right; i++){
        if(i != right) {
            int curr = (key[i].first >> (2*level)) & 0x3;
            if(prev == curr) continue;
            prev = curr;
        }

        // ここで [start, i) が1つの子ノード
        if(level > 0){
            int child_id = (key[start].first >> (2*level)) & 0x3;
            (*nid)++;
            node[id].children[child_id] = *nid;
            treeConstruction(node, particle, key, csum,
                            nid, level, start, i);
        }
        start = i;
    }

}


// pos, sizeを計算しながら再帰
void treeTrace(vector<TreeNode> &nodes, int node_id, double pos[2], double size) {
    TreeNode *node = &nodes[node_id];
    node->pos[0] = pos[0];
    node->pos[1] = pos[1];
    node->size = size;
    fprintf(stdout, "Node ID: %d\n", node->id);
    fprintf(stdout, "Center of Mass: (%.4f, %.4f)\n", node->cm[0], node->cm[1]);
    fprintf(stdout, "Node Pos: (%.4f, %.4f)\n", node->pos[0], node->pos[1]);
    fprintf(stdout, "Node Size: %.4f\n", node->size);
    fprintf(stdout, "Node Mass: %.4f\n\n", node->mass);
    double half = size * 0.5;
    for(int i=0; i<4; i++){
        if(node->children[i]==-1) continue;
        double child_pos[2];
        child_pos[0] = pos[0] + half * (i % 2);
        child_pos[1] = pos[1] + half * (i / 2);
        treeTrace(nodes, node->children[i], child_pos, half);
    }
}


void init_condition(vector<Particle> &particle, int n){
    for(int i=0; i<n; i++){
        particle[i].pos[0] = drand48();
        particle[i].pos[1] = drand48();
        particle[i].mass = 1.0/n;
    }
}

void init_csum(vector<Particle> &particle, vector<pair<int,int> > &key, vector<vector<double>> &csum, int n){
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
    int x0 = (int)(x * IMG_SIZE);
    int y0 = (int)(y * IMG_SIZE);
    int s = (int)(size * IMG_SIZE);
    // 上下
    for(int i = x0; i < x0 + s; i++) {
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
    for(int j = y0; j < y0 + s; j++) {
        if(x0 >= 0 && x0 < IMG_SIZE && j >= 0 && j < IMG_SIZE) {
            img[j][x0][0] = 0;
            img[j][x0][1] = 255;
            img[j][x0][2] = 0;
        }
        if(x0 + s - 1 >= 0 && x0 + s - 1 < IMG_SIZE && j >= 0 && j < IMG_SIZE) {
            img[j][x0 + s - 1][0] = 0;
            img[j][x0 + s - 1][1] = 255;
            img[j][x0 + s - 1][2] = 0;
        }
    }
}

void drawTree(unsigned char img[IMG_SIZE][IMG_SIZE][3], vector<TreeNode> &nodes, int node_id) {
    TreeNode *node = &nodes[node_id];
    drawBoundary(img, node->pos[0], node->pos[1], node->size);
    for(int i = 0; i < 4; i++) {
        if(node->children[i] != -1) {
            drawTree(img, nodes, node->children[i]);
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
    drawTree(img, node, root_id);
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
    if(n > MAX_PARTICLES){
        fprintf(stderr, "Error: the number of particles is too large");
        exit(1);
    }

    // levelの扱いがバラバラになってるので見直しが必要. 
    // levelの数の階層までしか掘り下げられない.
    // levelはツリーの深さを決めるパラメータ。大きいほど細かく分割されるが、計算量も使用メモリも増える。
    // levelに応じて、最下層は4^(level-1)分割されることになる。
    int nid = 0;
    int level = MAX_LEVEL; 
    int left = 0, right = n;

    vector<Particle> particle(n);
    vector<pair<int,int> > key(n);
    vector<TreeNode> nodes(NMAX);
    vector<vector<double>> csum(n+1, vector<double>(3, 0.0)); // 累積和配列(0: x*mass, 1: y*mass, 2: mass)
    
    init_condition(particle, n);
    p2key(particle, key, n, level);
    init_csum(particle, key, csum, n); 

    // ツリー構築
    treeConstruction(nodes, particle, key, csum, &nid, level, left, right);

    // pos/sizeを再帰的にセットしながらトレース
    double root_pos[2] = {0.0, 0.0};
    double root_size = 1.0;
    treeTrace(nodes, 0, root_pos, root_size);

    // 画像出力
    outputBitmap("output.ppm", nodes, 0, particle, n);
    return 0;
}
