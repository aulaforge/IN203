#include <cstdlib>
#include <string>
#include <iostream>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <fstream>
#include <ctime>
#include <iomanip>      // std::setw
#include <chrono>
#include <thread>
#include <mpi.h>
#include <unistd.h>
#include <stdint.h>

#include "parametres.hpp"
#include "galaxie.hpp"

;

//pour activer la méthode omp, "désannoter" le # pragma omp dans le fichier parametre.cpp

int main(int argc, char ** argv)
{

    char method[]="process"; //choix entre méthode en mémoire partagée et en mémoire distribuée (process==méthode MPI, thread==méthode thread)
    int MesureSpeedUp=0;  //calcul du speed up=variable à 1, sinon égale à 0
    if(strcmp(method,"thread")==0){
      std::cout<<"thread"<<std::endl; //affichage en début de programme pour s'assurer du bon choix de méthode, pause de 2 secondes pour laisser le temps de lire le texte
      sleep(2);
      char commentaire[4096];
      int width, height;
      SDL_Event event;
      SDL_Window   * window;

      parametres param;


      std::ifstream fich("parametre.txt");
      fich >> width;
      fich.getline(commentaire, 4096);
      fich >> height;
      fich.getline(commentaire, 4096);
      fich >> param.apparition_civ;
      fich.getline(commentaire, 4096);
      fich >> param.disparition;
      fich.getline(commentaire, 4096);
      fich >> param.expansion;
      fich.getline(commentaire, 4096);
      fich >> param.inhabitable;
      fich.getline(commentaire, 4096);
      fich.close();

      std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
      std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
      std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
      std::cout << "\t Chance expansion : " << param.expansion << std::endl;
      std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
      std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;
      std::srand(std::time(nullptr));

      SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);

      window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                width, height, SDL_WINDOW_SHOWN);

      galaxie g(width, height, param.apparition_civ);
      galaxie g_next(width, height);
      galaxie_renderer gr(window);

      int deltaT = (20*52840)/width;
      std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

      std::cout << std::endl;

      //gr.render(g); code de base
      unsigned long long temps = 0;

      std::chrono::time_point<std::chrono::system_clock> start, end1; //end2;
      int nbTour=0;
      float MoyenneTemps=0;
      while (1) {

          start = std::chrono::system_clock::now();
          std::thread t1(mise_a_jour,param, width, height, g.data(), g_next.data()); //lancement d'un thread parallèle pour mettre à jour les paramètres
          gr.render(g); //affichage du "tour" précédent pendant le calcul de la période suivante
          //end1 = std::chrono::system_clock::now();
          t1.join(); //mise à jour effective via méthode thread
          g_next.swap(g);
          //gr.render(g); code de base
          //end2 = std::chrono::system_clock::now();
          end1 = std::chrono::system_clock::now();
          std::chrono::duration<double> elaps1 = end1 - start;
        //std::chrono::duration<double> elaps2 = end2 - start;
          if (MesureSpeedUp==1){ //"cellule" de calcul du speed up
            nbTour++;
            MoyenneTemps+=elaps1.count()*1000;
            if (nbTour==1000){
              MoyenneTemps=MoyenneTemps/1000;
              std::cout<<std::setprecision(3)<<"Moyenne du temps par tour = "<<MoyenneTemps<<std::endl;
              return EXIT_SUCCESS;
            }
          }


          temps += deltaT;
          std::cout << "Temps passe : "
                  << std::setw(10) << temps << " années"
                  << std::fixed << std::setprecision(3)
                //<< "  " << "|  CPU(ms) : affichage " << elaps1.count()*1000
                  << "  " << "| CPU(ms) : affichage+calcul (recouvrement entrée/sortie) " << elaps1.count()*1000 //un seul affichage de temps car l'affichage et le calcul sont fait en parallèle
                  << "\r" << std::flush;
                  //_sleep(1000);
                  if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                    std::cout << std::endl << "The end" << std::endl;
                    break;
                  }
        }
        SDL_DestroyWindow(window);
        SDL_Quit();
        return EXIT_SUCCESS;
    }

    else if (strcmp(method,"process")==0){ //affichage en début de programme pour s'assurer du bon choix de méthode, pause de 2 secondes pour laisser le temps de lire le texte
      std::cout<<"process"<<std::endl;
      sleep(2);
      MPI_Init(&argc,&argv); //initialisation de MPI
      MPI_Status status;
      int rank, numtasks;
      MPI_Comm_size ( MPI_COMM_WORLD , & numtasks );

      if (numtasks<2){ //si un seul processus, exécution du programme de base
        char commentaire[4096];
        int width, height;
        SDL_Event event;
        SDL_Window   * window;

        parametres param;


        std::ifstream fich("parametre.txt");
        fich >> width;
        fich.getline(commentaire, 4096);
        fich >> height;
        fich.getline(commentaire, 4096);
        fich >> param.apparition_civ;
        fich.getline(commentaire, 4096);
        fich >> param.disparition;
        fich.getline(commentaire, 4096);
        fich >> param.expansion;
        fich.getline(commentaire, 4096);
        fich >> param.inhabitable;
        fich.getline(commentaire, 4096);
        fich.close();

        std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
        std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
        std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
        std::cout << "\t Chance expansion : " << param.expansion << std::endl;
        std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
        std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;
        std::srand(std::time(nullptr));

        SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);

        window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                  width, height, SDL_WINDOW_SHOWN);

        galaxie g(width, height, param.apparition_civ);
        galaxie g_next(width, height);
        galaxie_renderer gr(window);

        int deltaT = (20*52840)/width;
        std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

        std::cout << std::endl;

        gr.render(g);
        unsigned long long temps = 0;

        std::chrono::time_point<std::chrono::system_clock> start, end1, end2;
        float MoyenneTemps=0;
        int nbTour=0;
        while (1) {
            start = std::chrono::system_clock::now();
            mise_a_jour(param, width, height, g.data(), g_next.data());
            end1 = std::chrono::system_clock::now();
            g_next.swap(g);
            gr.render(g);
            end2 = std::chrono::system_clock::now();

            std::chrono::duration<double> elaps1 = end1 - start;
            std::chrono::duration<double> elaps2 = end2 - end1;
            std::chrono::duration<double> ttot =end2 - start;

            if (MesureSpeedUp==1){
              nbTour++;
              MoyenneTemps+=ttot.count()*1000;
              if (nbTour==1000){
                MoyenneTemps=MoyenneTemps/1000;
                std::cout<<std::setprecision(3)<<"Moyenne du temps par tour = "<<MoyenneTemps<<std::endl;
                MPI_Finalize();
                return EXIT_SUCCESS;
              }
            }


            temps += deltaT;
            std::cout << "Temps passe : "
                      << std::setw(10) << temps << " années"
                      << std::fixed << std::setprecision(3)
                      << "  " << "|  CPU(ms) : calcul " << elaps1.count()*1000
                      << "  " << "affichage " << elaps2.count()*1000
                      << "\r" << std::flush;
            //_sleep(1000);
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
              std::cout << std::endl << "The end" << std::endl;
              break;
            }
        }
        SDL_DestroyWindow(window);
        SDL_Quit();

        MPI_Finalize();
        return EXIT_SUCCESS;
      }



      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      char commentaire[4096];
      int width, height;
      SDL_Event event;
      SDL_Window   * window;

      parametres param;


      std::ifstream fich("parametre.txt");
      fich >> width;
      fich.getline(commentaire, 4096);
      fich >> height;
      fich.getline(commentaire, 4096);
      fich >> param.apparition_civ;
      fich.getline(commentaire, 4096);
      fich >> param.disparition;
      fich.getline(commentaire, 4096);
      fich >> param.expansion;
      fich.getline(commentaire, 4096);
      fich >> param.inhabitable;
      fich.getline(commentaire, 4096);
      fich.close();
      if(rank==0){
        std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
        std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
        std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
        std::cout << "\t Chance expansion : " << param.expansion << std::endl;
        std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
        std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;



      }
      SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);

      window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                width, height, SDL_WINDOW_SHOWN);
      std::srand(std::time(nullptr));
      galaxie g(width, height, param.apparition_civ);
      galaxie g_next(width, height);
      galaxie_renderer gr(window);

      int deltaT = (20*52840)/width;
      if(rank==0){
        std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

        std::cout << std::endl;
      }
      //gr.render(g); code de base
      unsigned long long temps = 0;


      std::chrono::time_point<std::chrono::system_clock> start, end1; //end2;
      int nbTour=0;
      float MoyenneTemps=0;
      if (numtasks==2){ //ici simple recouvrement entrée sortie
        while (1) {
            if (rank==1){ //reçoit les données du processus maitre (celui qui affiche), calcule l'état suivant et le renvoie au processus maitre
              MPI_Recv(g.data(), width*height, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              mise_a_jour(param, width, height, g.data(), g_next.data());
              MPI_Send(g_next.data(), width*height, MPI_CHAR, 0, 101, MPI_COMM_WORLD);
            }
            else if(rank==0){
              start = std::chrono::system_clock::now();
              MPI_Send(g.data(), width*height, MPI_CHAR, 1, 100, MPI_COMM_WORLD); //envoie les données au processus 1 pour qu'il calcule
              gr.render(g); //affichage du "tour" précédent pendant le calcul de la période suivante
              MPI_Recv(g_next.data(), width*height, MPI_CHAR, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //réception du prochain état, calculé par le processus 1
              g_next.swap(g); //mise à jour de l'état de la galaxie avec l'état calculé


              end1 = std::chrono::system_clock::now();
              std::chrono::duration<double> elaps1 = end1 - start;
            //std::chrono::duration<double> elaps2 = end2 - start;

              if (MesureSpeedUp==1){
                nbTour++;
                MoyenneTemps+=elaps1.count()*1000;
                if (nbTour==1000){
                  MoyenneTemps=MoyenneTemps/1000;
                  std::cout<<std::setprecision(3)<<"Moyenne du temps par tour = "<<MoyenneTemps<<std::endl;
                  MPI_Finalize();
                  return EXIT_SUCCESS;
                }
              }

              temps += deltaT;
              std::cout << "Temps passe : "
                        << std::setw(10) << temps << " années"
                        << std::fixed << std::setprecision(3)
                      //<< "  " << "|  CPU(ms) : affichage " << elaps1.count()*1000
                        << "  " << "| CPU(ms) : affichage+calcul (recouvrement entrée/sortie) " << elaps1.count()*1000 //un seul affichage de temps car recouvrement entrée sortie, seul le temps total est pertinent
                        << "\r" << std::flush;
                        //_sleep(1000);
                        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                          std::cout << std::endl << "The end" << std::endl;
                          break;
                        }
            }

        }
      }
      else { //cas où on a plus de 2 processus
        int portion=height/(numtasks-1); //taille en terme de "hauteur" des portions que vont calculer les processus dans la galaxie totale, le dernier calcule en plus les lignes restantes (division entière ici donc il manque des lignes potentiellement) (les lignes de cellules fantômes ne sont pas dans cette hauteur, elles sont rajoutées au "cas par cas", selon l'mplacement du processus)
        char reception[((width*(portion+1)+width*(portion+2)*(numtasks-3)+width*(portion+height%(numtasks-1)+1))+1)]=""; //tableau qui va recevoir successivement les "sous galaxies" qu'auront calculé les processus, cellules fantômes comprises

        while (1) {
            if(rank==0){
              MPI_Request request[2*(numtasks-1)];
              for(int k=0; k<2*(numtasks-1); k++){
                request[k]=MPI_REQUEST_NULL;
              }
              //MPI_Status status[2*(numtasks-1)];
              start = std::chrono::system_clock::now();
              int tag=100;

              //envoi de la partie à calculer aux différents processus
              //choix d'envoi/réception asynchrone pour éviter qu'un processus plus lent ne ralentisse le rythme et empêche d'autres processus prêts à envoyer leurs données de le faire

              MPI_Isend(g.data(), width*(portion+1), MPI_CHAR, 1, tag, MPI_COMM_WORLD, &request[0]); //processus 1;
              for(int i=2; i<numtasks-1; i++){ //processus entre 2 et le nombre de processus-2 (bornes comprises)
                tag++;;
                MPI_Isend(&g.data()[width*(portion*(i-1)-1)], width*(portion+2), MPI_CHAR, i, tag, MPI_COMM_WORLD, &request[i-1]);
              }
              tag++;
              MPI_Isend(&g.data()[width*(portion*(numtasks-2)-1)], width*(portion+height%(numtasks-1)+1), MPI_CHAR, numtasks-1, tag, MPI_COMM_WORLD,&request[numtasks-2]); //dernier processus
              gr.render(g); //affichage du "tour" précédent pendant le calcul de la période suivante

              //réception du calcul effectué par les différents processus

              int indice=0;
              MPI_Irecv(reception, width*(portion+1), MPI_CHAR, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &request[numtasks-1]); //processus 1
              indice+=width*(portion+1);
              for(int i=2; i<numtasks-1; i++){ //processus entre 2 et le nombre de processus-2 (bornes comprises)
                MPI_Irecv(&reception[indice], width*(portion+2), MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &request[numtasks+i-2]);
                indice+=width*(portion+2);
              }
              MPI_Irecv(&reception[indice], width*(portion+height%(numtasks-1)+1), MPI_CHAR, numtasks-1, MPI_ANY_TAG, MPI_COMM_WORLD, &request[numtasks-3]); //dernier processus


              MPI_Waitall(2*(numtasks-1), request, NULL); //barrière car on a ensuite besoin de l'ensemble des données

              //comparaison et mise à jour des planètes "réelles" et des planètes "fantômes"

              indice=width*(portion+1);
              int pas=width*(portion+2);
              char planete_reelle, planete_fantome;


//phase de "traitement", on compare les valeurs entre les lignes fantômes et les lignes réelles correspondantes (les 1/0/-1 correspondent à habitée, habitable et inhabitable, j'avais modifié pour un test, et comme ça ne changeait pas mes résultats j'ai laissé comme ça pour éviter de refaire le travail inverse (ctrl z était impossible))

              //traitement du calcul du processus 1
              for(int j=0; j<width; j++){
                memcpy(&planete_reelle,&reception[width*(portion-1)+j], sizeof(char));
                memcpy(&planete_fantome,&reception[width*(portion+1)+j], sizeof(char));
                if((planete_reelle==0&&planete_fantome==1)||(planete_reelle==1&&planete_fantome==0)) {reception[width*(portion-1)+j]='1';}
                else if((planete_reelle==0&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==0)) {reception[width*(portion-1)+j]=-1;}
                else if((planete_reelle==1&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==1)) {reception[width*(portion-1)+j]=-1;}
              }
              memcpy(g_next.data(), reception, width*portion*sizeof(char)); //copie de la partie de la réception


              //traitement du calcul des processus entre 2 et le nombre de processus-2 (bornes comprises)
              for(int i=2; i<numtasks-1; i++){
                for(int j=0; j<width; j++){
                  memcpy(&planete_reelle,&reception[indice+width+j], sizeof(char));
                  memcpy(&planete_fantome,&reception[indice-width+j], sizeof(char));
                  if((planete_reelle==0&&planete_fantome==1)||(planete_reelle==1&&planete_fantome==0)) {reception[indice+width+j]='1';}
                  else if((planete_reelle==0&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==0)) {reception[indice+width+j]=-1;}
                  else if((planete_reelle==1&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==1)) {reception[indice+width+j]=-1;}
                }
                for(int j=0; j<width; j++){
                  memcpy(&planete_reelle,&reception[indice+pas-width+j],sizeof(char));
                  memcpy(&planete_fantome,&reception[indice+pas+width+j],sizeof(char));
                  if((planete_reelle==0&&planete_fantome==1)||(planete_reelle==1&&planete_fantome==0)) {reception[indice+pas-width+j]='1';}
                  else if((planete_reelle==0&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==0)) {reception[indice+pas-width+j]=-1;}
                  else if((planete_reelle==1&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==1)) {reception[indice+pas-width+j]=-1;}
                }
                memcpy(&g_next.data()[(i-1)*portion*width], &reception[indice+width], width*portion*sizeof(char));
                indice+=pas;
              }
              //traitement du calcul du processus nombre de processus-1
              for(int j=0; j<width; j++){
                memcpy(&planete_reelle,&reception[indice+width+j],sizeof(char));
                memcpy(&planete_fantome,&reception[indice-width+j],sizeof(char));
                if((planete_reelle==0&&planete_fantome==1)||(planete_reelle==1&&planete_fantome==0)) {reception[indice+width+j]='1';}
                else if((planete_reelle==0&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==0)) {reception[indice+width+j]=-1;}
                else if((planete_reelle==1&&planete_fantome==-1)||(planete_reelle==-1&&planete_fantome==1)) {reception[indice+width+j]=-1;}
                //std::cout<<planete_reelle<<std::endl;
                //sleep(4);
              }
              memcpy(&g_next.data()[(numtasks-2)*portion*width], &reception[indice+width], width*(portion+height%(numtasks-1))*sizeof(char));

              //std::cout<<g_next.data()[(numtasks-2)*portion*width]<<"            ///////                   "<<reception[indice+width]<<"  test: "<<habitable<<std::endl;
              //sleep(4);
              g_next.swap(g);


              end1 = std::chrono::system_clock::now();
              std::chrono::duration<double> elaps1 = end1 - start;
            //std::chrono::duration<double> elaps2 = end2 - start;
              if (MesureSpeedUp==1){
                nbTour++;
                MoyenneTemps+=elaps1.count()*1000;
                if (nbTour==1000){
                  MoyenneTemps=MoyenneTemps/1000;
                  std::cout<<std::setprecision(3)<<"Moyenne du temps par tour = "<<MoyenneTemps<<std::endl;
                  MPI_Finalize();
                  return EXIT_SUCCESS;
                }
              }

              temps += deltaT;
              std::cout << "Temps passe : "
                        << std::setw(10) << temps << " années"
                        << std::fixed << std::setprecision(3)
                      //<< "  " << "|  CPU(ms) : affichage " << elaps1.count()*1000
                        << "  " << "| CPU(ms) : affichage+calcul (recouvrement entrée/sortie) " << elaps1.count()*1000
                        << "\r" << std::flush;
                        //_sleep(1000);
                        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                          std::cout << std::endl << "The end" << std::endl;
                          break;
                        }
            }
            int tag=100+2*numtasks;
            if(rank==1){
              char ancienne_data[width*(portion+1)]; //buffer de récéption
              char nouv_data[width*(portion+1)]; //buffer qui accueillera les données calculées et qui sera envoyé
              MPI_Recv(ancienne_data, width*(portion+1), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              mise_a_jour(param, width, portion+1, ancienne_data, nouv_data); //calcul de l'état suivant de la partie associée au processus
              MPI_Send(nouv_data, width*(portion+1), MPI_CHAR, 0, tag+rank, MPI_COMM_WORLD);
            }
            for(int test=2; test<numtasks-1; test++){
              if (rank==test){
                char ancienne_data[width*(portion+2)]; //buffer de récéption
                char nouv_data[width*(portion+2)]; //buffer qui accueillera les données calculées et qui sera envoyé
                MPI_Recv(ancienne_data, width*(portion+2), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                mise_a_jour(param, width, portion+2, ancienne_data, nouv_data); //calcul de l'état suivant de la partie associée au processus
                MPI_Send(&nouv_data, width*(portion+2), MPI_CHAR, 0, tag+rank, MPI_COMM_WORLD);
              }
            }
            if(rank==numtasks-1){
              char ancienne_data[width*(portion+height%(numtasks-1)+1)]; //buffer de récéption
              char nouv_data[width*(portion+height%(numtasks-1)+1)]; //buffer qui accueillera les données calculées et qui sera envoyé
              MPI_Recv(ancienne_data, width*(portion+height%(numtasks-1)+1), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              mise_a_jour(param, width, portion+height%(numtasks-1)+1, ancienne_data, nouv_data); //calcul de l'état suivant de la partie associée au processus
              MPI_Send(&nouv_data, width*(portion+height%(numtasks-1)+1), MPI_CHAR, 0, tag+rank, MPI_COMM_WORLD);
            }
        }
      }
      if(rank==0){
        SDL_DestroyWindow(window);
        SDL_Quit();
      }
      MPI_Finalize();
      return EXIT_SUCCESS;
    }
}
