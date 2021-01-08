/************************************************************************/
/* 重要變數與參數說明                                                   */
/*                                                                      */
/* 1.particle_cnt: 粒子數, 通常於-40,複雜設-200                         */
/*                                                                      */
/* 2.粒子長度: 即解空間維度,或使用變數個數,由於此問題之適應函式只用到   */
/*             x 變數, 故粒子長度為，將再補程式碼做為考慮粒子長度為n    */
/*             之情形                                                   */
/*                                                                      */
/* 3.max_pos, min_pos: 粒子範圍,即解空間中,各維度(變數)之範圍限制,      */
/*             普遍性應考慮n 維之情形                                   */
/*                                                                      */
/* 4.max_v : 即最大速度限制, 通常設定成粒子之寬度,即解空間範圍,普遍性考 */
/*           慮n 維情形,即每維之解空間範圍                              */
/* 5.c1,c2 : 學習常數, c1,c2多設, 一般c1=c2=[0,4]                       */
/*                                                                      */
/* 6.w     : 慣性權重, 此常數為後來由Yuhui Shi and Russell Eberhart 提  */
/*           出,多設在[0,1.5] 之間, 設1.0 為保守值                      */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------
// fitness value
// best value on [-100, 100] , aboue fit(-57.469) = 390245.791738

double fit(double x, double y)
{
     
     // x**3 - 0.8x**2 - 1000x + 8000
     // return fabs(8000.0 + x*(-10000.0+x*(-0.8+x)));

     return pow((x - 20), 2) + pow((y - 20), 2) + 1;
}


//////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------
// * 使用pso 求x**3 - 0.8x**2 - 10000x + 8000 於[-100,100] 間之最大值    *
// -----------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////

/* 定義結構體particle */
typedef struct tag_particle{
     double x;  /* 目前位置, 即x value    */
     double y;
     double velocity_x;  /* 目前粒子速度           */
     double velocity_y;
     double fitness ;  /* 適應函式值             */
     double pbest_pos_x; /* particle 目前最好位置  */
     double pbest_pos_y;
     double pbest_fit; /* particle 目前最佳適應值*/
}particle;

// -----------------------------------------------------------------------
// * 全域變數宣告, 較好佈局為拆pos.h / pos.c , 放在pos.c 宣告成static    *
// * 再做set_param function 為入口, 設定所有static param variable        *
// * 或直接用#define 方式, 也需get_gbest 做入口取得全域最佳值            *
// -----------------------------------------------------------------------

double w, c1, c2;                    /* 相關權重參數   */
double max_v;                        /* 最大速度限制   */
double max_pos, min_pos;             /* 最大,小位置限制*/
int particle_cnt;               /* 粒子數量       */
particle gbest;                      /* 全域最佳值     */

// -----------------------------------------------------------------------
// * pso 相關函式宣告                                                    *
// -----------------------------------------------------------------------

#define RND() ((double)rand()/RAND_MAX) /* 產生[0,1] 亂數         */
particle* AllocateParticle();           /* 配置particle_cnt 個粒子*/
void ParticleInit(particle *p);         /* 粒子初始化             */
void ParticleMove(particle *p);         /* 開始移動               */
void ParticleRelease(particle* p);      /* 釋放粒子記憶體         */
void ParticleDisplay(particle* p);      /* 顯示所有粒子資訊       */

//////////////////////////////////////////////////////////////////////////

int main()
{
    /* 變數宣告*/
    unsigned i;
    unsigned max_itera = 100;           /* max_itera : 最大演化代數*/
    particle* p;                         /* p         : 粒子群      */

    /* 設定參數*/
    min_pos = 0.f , max_pos = 50.f;  /* 位置限制, 即解空間限制   */
    w = 0.2, c1 = 0.1, c2 = 0.1;          /* 慣性權重與加速常數設定   */
    particle_cnt = 100;                  /* 設粒子個數               */
    max_v = (max_pos-min_pos) * 1.0;    /* 設最大速限               */

    /* 開始進行*/
    float start = clock();
    p = AllocateParticle();      // 配置記憶體
    ParticleInit(p);             // 粒子初始化
    for(i = 0; i < max_itera; i++)   // 進行迭代
        ParticleMove(p);       // 粒子移動
    ParticleDisplay(p);          // 顯示最後結果
    float finish = clock();
    double duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Executing time: %f\n", duration);
    ParticleRelease(p);          // 釋放記憶體

    return 0;
}

//////////////////////////////////////////////////////////////////////////

// 配置particle_cnt 個粒子
particle* AllocateParticle()
{
     return (particle*)malloc(sizeof(particle)*particle_cnt);
}

// 粒子初始化
void ParticleInit(particle *p)
{
     unsigned i;
     const double pos_range = max_pos - min_pos; // 解寬度
     srand((unsigned)time(NULL));

     // 以下程式碼效率不佳, 但較易懂一點
     for(i = 0; i < particle_cnt; i++) {
          // 隨機取得粒子位置, 並設為該粒子目前最佳適應值
          p[i].pbest_pos_x = p[i].x = RND() * pos_range + min_pos; 
          p[i].pbest_pos_y = p[i].y = RND() * pos_range + min_pos; 
          // 隨機取得粒子速度
          p[i].velocity_x = RND() * max_v;
          p[i].velocity_y = RND() * max_v;
          // 計算該粒子適應值, 並設為該粒子目前最佳適應值
          p[i].pbest_fit = p[i].fitness = fit(p[i].x, p[i].y);

          // 全域最佳設定
          if(i==0 || p[i].pbest_fit > gbest.fitness) 
               memcpy((void*)&gbest, (void*)&p[i], sizeof(particle));
     }
}

// 開始移動
void ParticleMove(particle *p)
{
     unsigned i;
     double v_x, v_y, pos_x, pos_y;     // 暫存每個粒子之速度, 位置用
     double ppos_x, ppos_y, gpos_x, gpos_y; // 暫存區域及全域最佳位置用
     gpos_x = gbest.x;
     gpos_y = gbest.y;

     // 更新速度與位置
     for(i = 0; i < particle_cnt; i++){
           v_x = p[i].velocity_x; // 粒子目前速度
           v_y = p[i].velocity_y;
           pos_x = p[i].x; // 粒子目前位置
           pos_y = p[i].y;
           ppos_x = p[i].pbest_pos_x; // 粒子目前曾到到最好位置
           ppos_y = p[i].pbest_pos_y;
           
           v_x = w * v_x + c1*RND()*(ppos_x - pos_x) + c2*RND()*(gpos_x - pos_x); // 更新速度
           v_y = w * v_y + c1*RND()*(ppos_y - pos_y) + c2*RND()*(gpos_y - pos_y);
          
          
          pos_x = pos_x + v_x;               // 更新位置
          pos_y = pos_y + v_y;

          p[i].velocity_x = v_x;        // 更新粒子速度   
          p[i].velocity_y = v_y;   
          p[i].x = pos_x;       // 更新粒子位置
          p[i].y = pos_y;
          p[i].fitness = fit(pos_x, pos_y); // 更新粒子適應值

          // 更新該粒子目前找過之最佳值
          if(p[i].fitness < p[i].pbest_fit) {
               p[i].pbest_fit = p[i].fitness;
               p[i].pbest_pos_x = p[i].x;
               p[i].pbest_pos_y = p[i].y;
          }

          // 更新全域最佳值
          if(p[i].fitness < gbest.fitness) 
               memcpy((void*)&gbest, (void*)&p[i], sizeof(particle));
     }
}

// 釋放粒子記憶體
void ParticleRelease(particle* p)
{
     free(p);
}

// 顯示所有粒子資訊
void ParticleDisplay(particle* p)
{
     unsigned i;

     /* 若想看完整的粒子資料，可把下面三行註解拿掉，這裡只顯示最佳解。*/
     // for(i=0; i<particle_cnt; i++)
     //      printf("#%d : %lf , %lf, %lf\n", i+1, p[i].x, p[i].y, p[i].fitness);
     // puts("------------------------------\n");
     printf("best : %10.6lf , %10.6lf, %lf\n", gbest.x, gbest.y, gbest.fitness);     
}