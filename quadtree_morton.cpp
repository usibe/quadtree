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
    int id;
    double cm[2];
    double pos[2];
    double mass;
    double size;
    TreeNode *children[4];
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
                     int *nid, int level, int left, int right){
    
    level -= 1;
    double n = right - left;
    int id = *nid;
    for(int i=left; i<right; i++){
        int j = key[i].second;
        node[id].cm[0] += particle[j].pos[0];
        node[id].cm[1] += particle[j].pos[1];
        node[id].mass += particle[j].mass;
    }
    node[id].cm[0]/=n;
    node[id].cm[1]/=n;
    node[id].id = id; //いるかこれ？

    int npart[4] = {0, 0, 0, 0};
    for(int i=left; i<right; i++){
        int ii = (key[i].first >> (2*level)) & 0x3;
        npart[ii]++;
    }

    right = left;
    for(int i=0; i<4; i++){
        if(npart[i] == 0 || level<=0){
            node[id].children[i] = nullptr;
            continue;
        }
        *nid+=1;
        left = right;
        right += npart[i];
        node[id].children[i] = &node[*nid]; //childrenに入るのは配列番号だけで良いのでは？
        treeConstruction(node, particle, key, nid,
                         level, left, right);
    }

}

/*
//ノードを最初からセル分割しておいて、ツリー構築を行わないパターン
void treeConstruction(vector<TreeNode> &node, vector<Particle> &particle,
                     vector<int> &index, int level, int n){

    int nid = 0;
    //各粒子のmortonkeyからノードを取得。親の方から順番にposとmassを足していく。最後に平均とって重心出す。
    for(int i=0; i<n; i++){
        
        int x, y, scale=1<<level;
        x = scale * particle[i].pos[0];
        y = scale * particle[i].pos[1];
        int key = get_key(x, y, level);

        for(int j=level-1; j>=0; j--){
            //mortonkeyから頑張って辿る
            //子セルがnullじゃなければセルを作成。そうじゃなければ次の階層へスキップ。
            //タイムステップごとの更新が面倒になるな...どうしようか...
            int id = (key >> (2*j)) & 0x3;
            node[nid].cm[0] += particle[i].pos[0];
            node[nid].cm[1] += particle[i].pos[1];
            node[nid].mass += particle[i].mass;
            if() node[nid].children[id] = id; //子ノードのアドレス登録　idはmortonkeyから取ってくる必要がある。
            //nodeのchildren配列はnull初期化しないといけないな...
        }
    }
}
*/


// pos, sizeを計算しながら再帰
void treeTrace(TreeNode *node, double pos[2], double size) {
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
        if(node->children[i]==nullptr) continue;
        double child_pos[2];
        child_pos[0] = pos[0] + half * (i % 2);
        child_pos[1] = pos[1] + half * (i / 2);
        treeTrace(node->children[i], child_pos, half);
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

void drawTree(unsigned char img[IMG_SIZE][IMG_SIZE][3], TreeNode *node) {
    drawBoundary(img, node->pos[0], node->pos[1], node->size);
    for(int i = 0; i < 4; i++) {
        if(node->children[i]) {
            drawTree(img, node->children[i]);
        }
    }
}

void outputBitmap(const char *filename, TreeNode *root, vector<Particle> &particle, int n) {
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
    drawTree(img, root);
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
    vector<TreeNode> node(NMAX);
    
    init_condition(particle, n);
    p2key(particle, key, n, level); 

    // ツリー構築
    treeConstruction(node, particle, key, &nid, level, left, right);

    // pos/sizeを再帰的にセットしながらトレース
    double root_pos[2] = {0.0, 0.0};
    double root_size = 1.0;
    treeTrace(&node[0], root_pos, root_size);

    // 画像出力
    outputBitmap("output.ppm", &node[0], particle, n);
    return 0;
}
