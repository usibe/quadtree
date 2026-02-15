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


void treeConstruction(vector<TreeNode> &node, int *pcid, vector<Particle> &particle,
                     vector<int> &index, vector<int> &index_tmp, int n, int fst){
    int cid = *pcid;
    double size = node[cid].size*0.5;
    int npart[4] = {0,0,0,0};
    int child_fst[4] = {0,0,0,0};
    int child_lst[4] = {0,0,0,0};
    int ntot = 0;
    double cmtot[2] = {0,0};
    
    node[cid].id = cid;

    if(n == 1){
        node[cid].cm[0] = particle[fst].pos[0];
        node[cid].cm[1] = particle[fst].pos[1];
        node[cid].mass = particle[fst].mass;
        return;
    }

    for(int i=fst; i<fst+n; i++){
        int id = index[i];
        int i0 = (particle[id].pos[0] - node[cid].pos[0])/size;
        int i1 = (particle[id].pos[1] - node[cid].pos[1])/size;
        int cid = i0 + 2*i1;
        npart[cid]++;
    }

    for(int i=0; i<4; i++){
        child_fst[i] = ntot + fst;
        ntot += npart[i];
    }

    for(int i=fst; i<fst+n; i++){
        int id = index[i];
        int i0 = (particle[id].pos[0] - node[cid].pos[0])/size;
        int i1 = (particle[id].pos[1] - node[cid].pos[1])/size;
        int cid = i0 + 2*i1;
        int j = child_fst[cid] + child_lst[cid];
        child_lst[cid]++;
        index_tmp[j] = id;
    }

    for(int i=0; i<n; i++) index[fst+i] = index_tmp[fst+i];

    for(int i=0; i<4; i++){
        if(npart[i] < 1){
            node[cid].children[i] = nullptr;
            continue;
        }

        *pcid+=1;
        int nid = *pcid;
        // int nid = 4*cid+i+1;
        if(*pcid >= NMAX){
            fprintf(stderr, "Error: there are too many nodes");
            exit(1);
        }
        
        node[cid].children[i] = &node[nid];
        node[nid].pos[0] = node[cid].pos[0] + size*(i%2);
        node[nid].pos[1] = node[cid].pos[1] + size*(int)(i/2);
        node[nid].cm[0] = 0.0;
        node[nid].cm[1] = 0.0;
        node[nid].mass = 0.0;
        node[nid].size = size;
        
        treeConstruction(node, pcid, particle, index, index_tmp, npart[i], child_fst[i]);
        node[cid].mass += node[nid].mass;
        node[cid].cm[0] += node[nid].cm[0]*node[nid].mass;
        node[cid].cm[1] += node[nid].cm[1]*node[nid].mass;
        cmtot[0] += node[nid].cm[0];
        cmtot[1] += node[nid].cm[1];
    }
    node[cid].cm[0] /= cmtot[0];
    node[cid].cm[1] /= cmtot[1];
}


void treeTrace(TreeNode *node){
    fprintf(stdout, "Node ID: %d\n", node->id);
    fprintf(stdout, "Center of Mass: (%.4f, %.4f)\n", node->cm[0], node->cm[1]);
    fprintf(stdout, "Node Pos: (%.4f, %.4f)\n", node->pos[0], node->pos[1]);
    fprintf(stdout, "Node Mass: %.4f\n\n", node->mass);
    
    for(int i=0; i<4; i++){
        if(node->children[i]==nullptr) continue;
        treeTrace(node->children[i]);
    }
}


void init_condition(vector<Particle> &particle, vector<int> &index, int n){
    for(int i=0; i<n; i++){
        particle[i].pos[0] = drand48();
        particle[i].pos[1] = drand48();
        particle[i].mass = 1.0/n;
        index[i] = i;
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

    vector<Particle> particle(n);
    vector<int> index(n), index_tmp(n);
    init_condition(particle, index, n);
    vector<TreeNode> node(NMAX); 
    
    node[0].pos[0] = 0.0;
    node[0].pos[1] = 0.0;
    node[0].size = 1.0;

    int cid = 0;
    treeConstruction(node, &cid, particle, index, index_tmp, n, 0);
    treeTrace(&node[0]);
    
    // 画像出力
    outputBitmap("output.ppm", &node[0], particle, n);
    return 0;
}