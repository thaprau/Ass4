
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <pthread.h>

typedef struct QuadTree {
        struct QuadTree *child1;
        struct QuadTree *child2;
        struct QuadTree *child3;
        struct QuadTree *child4;
        double mass;
        double C_O_x;
        double C_O_y;
        double posx;
        double posy;
        double side;
        struct particle ** particles;
} QT;

// Global variables
const double epsilon = 0.001;


// Struct desrcibing a particle
typedef struct particle {
    double x_force;
    double y_force;
    double x;
    double y;
    double mass;
    double x_vel;
    double y_vel;
    double brightness;
} part;

void force_x(part *p1, part *p2, double r, double G) {
 

    double force = - G * (p1->mass) * (p2->mass) * ((p1->x) - (p2->x)) /(
    (r + epsilon)*(r + epsilon)*(r + epsilon));

    p1->x_force = p1->x_force + force;
    p2->x_force = p2->x_force - force;

}

void force_y(part *p1, part *p2, double r, double G) {

    double force = - G * (p1->mass) * (p2->mass) * ((p1->y) - (p2->y)) / (
    (r + epsilon)*(r + epsilon)*(r + epsilon));

    p1->y_force = p1->y_force + force;
    p2->y_force = p2->y_force - force;
}

void vel_update(part** particles, int N, double t) {

    for(int i = 0; i < N; i++) {
        particles[i]->x_vel = particles[i]->x_vel + ((particles[i]->x_force)/(particles[i]->mass)) * t;
        particles[i]->y_vel = particles[i]->y_vel + ((particles[i]->y_force)/(particles[i]->mass)) * t;    
    }
}

void pos_update(part** particles, int N, double t) {

        for(int i = 0; i < N; i++) {
            particles[i]->x = particles[i]->x + (particles[i]->x_vel) * t;
            particles[i]->y = particles[i]->y + (particles[i]->y_vel) * t;
            particles[i]->y_force = 0;
            particles[i]->x_force = 0;

    }

}

void create_tree(QT *node, int N) {
        QT * child = (QT *)malloc(4*sizeof(QT));
        node->child1 = &child[0];
        node->child2 = &child[1];
        node->child3 = &child[2];
        node->child4 = &child[3];
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        double new_side = node->side/2;
        double posx = node->posx;
        double posy = node->posy;
        node->child1->side = new_side;
        node->child1->posx = posx;
        node->child1->posy = posy;
        node->child2->side = new_side;
        node->child2->posx = posx + new_side;
        node->child2->posy = posy;
        node->child3->side = new_side;
        node->child3->posx = posx;
        node->child3->posy = posy + new_side;
        node->child4->side = new_side;
        node->child4->posx = posx + new_side;
        node->child4->posy = posy + new_side;
        
        part **part1 = (part**)malloc(N*sizeof(part*));
        part **part2 = (part**)malloc(N*sizeof(part*));
        part **part3 = (part**)malloc(N*sizeof(part*));
        part **part4 = (part**)malloc(N*sizeof(part*));



        for(int i=0; i<N; i++) {
            if(node->particles[i]->y <= posy + new_side) {
                if(node->particles[i]->x <= posx + new_side) {
                    part1[count1] = node->particles[i];
                    count1++;
                } else {
                    part2[count2] = node->particles[i];
                    count2++;
                }
            } else {
                if(node->particles[i]->x <= posx + new_side) {
                    part3[count3] = node->particles[i];
                    count3++;    
                } else {
                    part4[count4] = node->particles[i];
                    count4++;
                }
            }
        }
        node->child1->particles = part1;
        node->child2->particles = part2;
        node->child3->particles = part3;
        node->child4->particles = part4;

        
        if(count1 > 1) {
            create_tree(node->child1, count1);
        }
        if(count2 > 1) {
            create_tree(node->child2, count2);
        }
        if(count3 > 1) {
            create_tree(node->child3, count3);
        }
        if(count4 > 1) {
            create_tree(node->child4, count4);
        }
}

void centerOFMass(QT * node, int N){
    double mass = 0;
    double x = 0;
    double y = 0;
    for(int i=0; i<N; i++) {
        mass += node->particles[i]->mass;
        x += node->particles[i]->mass*node->particles[i]->x;
        y += node->particles[i]->mass*node->particles[i]->y;
    }
    node->C_O_x = x/(mass);
    node->C_O_y = y/(mass);
    node->mass = mass;
    
}





int main(int argc, char* argv[]) {

    if (argc != 6) {
        printf("Wrong number of inputs");
        return -1;
    }
    // Read user input
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const double theta_max = atof(argv[5]);
    const int graphics = atoi(argv[6]);


    // Creates array to store particles
    struct particle **array = (part**) malloc(N*sizeof(part*));
    for(int i = 0; i < N; i++)
    {
        array[i] = (part*) malloc(sizeof(part));
    }

    // Read file
    
    FILE * fp = fopen(filename, "r");

    for(int i = 0; i < N; i++)
    {

        fread(&(array[i]->x), sizeof(double), 1, fp);
        fread(&(array[i]->y), sizeof(double), 1, fp);
        fread(&(array[i]->mass), sizeof(double), 1, fp);
        fread(&(array[i]->x_vel), sizeof(double), 1, fp);
        fread(&(array[i]->y_vel), sizeof(double), 1, fp);
        fread(&(array[i]->brightness), sizeof(double), 1, fp);

    }


    fclose(fp);





    // Prepare for the loop
    int t = 0;
    const double G = (double) 100/N;
    const float W = 1; 
    const float L = 1;

    // Loop with graphics
    if(graphics==1) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0,1);

        while(t <= nsteps){
            ClearScreen();
            printf("Loop nr %d \n", t);

            for(int i = 0; i < N; i++) 
            {
                array[i]->y_force = 0;
                array[i]->x_force = 0;
            }

            for(int i = 0; i < N; i++)  
            {
                
                DrawCircle(array[i]->x, array[i]->y, W, L, 0.01, 0);

                for(int j = i+1; j < N; j++)
                {
                    double r = sqrt(pow((array[i]->x)-(array[j]->x), 2) + pow((array[i]->y)-(array[j]->y), 2));
                    force_x(array[i], array[j], r, G );
                    force_y(array[i], array[j], r, G);

                }
            }


            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t);
            
            Refresh();
            usleep(3000);

            t+=1;
        }
        FlushDisplay();
        CloseDisplay();
    }

    // Loop without graphics
    else {


        while(t < nsteps){}
            //printf("Loop nr --- %d\n", t);
            for(int i = 0; i < N; i++)  
            {
                for(int j = i+1; j < N; j++)
                {
                    //printf("particle nr %d x : %lf\n", i, (array[i]->x-array[j]->x)*(array[i]->x-array[j]->x));
                    double r = sqrt(
                        (array[i]->x-array[j]->x)*(array[i]->x-array[j]->x) + (array[i]->y-array[j]->y)*(array[i]->y-array[j]->y)
                        );
                    //printf("r : %lf \n", r);
                    force_x(array[i], array[j], r, G);
                    force_y(array[i], array[j], r, G);
                }
            //printf("particle nr %d force y : %lf\n", i, array[i]->y_force);
            }
            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t);
            t +=1;
        }                                
    }

    // Writing to new file
    FILE * pw = fopen("result.gal", "w");

    for(int i = 0; i < N; i++)
    {
        fwrite(&(array[i]->x), sizeof(double), 1, pw);
        fwrite(&(array[i]->y), sizeof(double), 1, pw);
        fwrite(&(array[i]->mass), sizeof(double), 1, pw);
        fwrite(&(array[i]->x_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->y_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->brightness), sizeof(double), 1, pw);

    }
    
    fclose(pw);

    for(int i = 0; i < N; i++)
    {
        free(array[i]);
    }
    free(array);
return 0;
}