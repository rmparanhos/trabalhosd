#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

// MAX char table (ASCII)
#define MAX 256

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, int n, char *substr, int m) {

	int d[MAX];
	int i, j, k;

	// pre-processing
	for (j = 0; j < MAX; j++)
		d[j] = m + 1;
	for (j = 0; j < m; j++)
		d[(int) substr[j]] = m - j;

	// searching
	i = m - 1;
	while (i < n) {
		k = i;
		j = m - 1;
		while ((j >= 0) && (string[k] == substr[j])) {
			j--;
			k--;
		}
		if (j < 0)
			return k + 1;
		i = i + d[(int) string[i + 1]];
	}

	return -1;
}

FILE *fdatabase, *fquery, *fout, *foutparalelo;

void openfiles() {

	fdatabase = fopen("dna.in", "r+");
	if (fdatabase == NULL) {
		perror("dna.in");
		exit(EXIT_FAILURE);
	}

	fquery = fopen("query.in", "r");
	if (fquery == NULL) {
		perror("query.in");
		exit(EXIT_FAILURE);
	}

	fout = fopen("dna.out", "w");
	if (fout == NULL) {
		perror("fout");
		exit(EXIT_FAILURE);
	}

	foutparalelo = fopen("dnaparalelo.out", "w");
	if (foutparalelo == NULL) {
		perror("foutparalelo");
		exit(EXIT_FAILURE);
	}

}

void closefiles() {
	fflush(fdatabase);
	fclose(fdatabase);

	fflush(fquery);
	fclose(fquery);

	fflush(fout);
	fclose(fout);

	fflush(foutparalelo);
	fclose(foutparalelo);
}

void remove_eol(char *line) {
	int i = strlen(line) - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

void junta_saida(int np){
	int i;
	rewind(fquery);
	rewind(fdatabase);
	char desc_query[30];
	char desc_dna[50];
	int resultados[np];
	while(fgets(desc_query,30, fquery) != NULL){
		printf("\n%s\n",desc_query);
		while(fgets(desc_dna,50,fdatabase) != NULL){
			printf("\n%s\n",desc_dna);
			for(i=0;i<np;i++){
				char saida_arq[30];
				sprintf(saida_arq,"dnaparalelo%d.out",i);
				FILE *leitura_saida = fopen(saida_arq,"r");
				int resultado;
				fscanf(leitura_saida,"%d",&resultado);
				
				resultados[i] = resultado;
			}
			int menor = resultados[0];
			for(i=1;i<np;i++){
				if(menor>resultados[i]){
					if(resultados[i]!=-1){
						menor = resultados[i];
					}
				}
			}
			printf("\n%d\n",menor);
		}
		fgets(desc_query,30, fquery);
	}
}

void slice_str(const char * str, char * buffer, size_t start, size_t end)
{
    size_t j = 0;
    for ( size_t i = start; i <= end; ++i ) {
        buffer[j++] = str[i];
    }
    buffer[j] = 0;
}

int divide(char *string, int tam_string, char *substr, int tam_substring, char *desc_query, char *desc_dna) {
	int result;
	//int tam_string_procs_ini = (tam_string/nprocs) + tam_substring - 1;
	//int tam_string_procs_meio = (tam_string/nprocs) + tam_substring - 1 + tam_substring - 1;
	//int tam_string_procs_fim = (tam_string/nprocs) + tam_substring - 1;
	//printf("\n%d\n",tam_string_procs_ini);
	//printf("%d\n",tam_string_procs_meio);
	//printf("%d\n",tam_string_procs_fim);
	int meu_rank, np, origem, destino, tag = 0;
    char msg[100];
    MPI_Status status;

    
    MPI_Comm_rank(MPI_COMM_WORLD, &meu_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np); // número de processadores

	char *saida[20];
	sprintf(saida,"dnaparalelo%d.out",meu_rank);
	foutparalelo = fopen(saida, "a");

	int margem = tam_substring-1;
	int i;
	int pos=0;
	if(meu_rank == 0){
		int msg_ini;
		int msg_fim;
		//printf("\nstring\n");
		//printf(string);
		//printf("\n");
	

		for(i=0;i<np;i++){
			
			if(i==0){
				//printf("Ini %d fim %d\n",pos,tam_string/nprocs+margem);
				
				pos+=tam_string/np+margem;
				
			}
			else if(i<=(np-2)){
				//printf("Ini %d fim %d\n",pos-margem,pos+tam_string/nprocs);
				msg_ini = pos-margem;
				msg_fim = pos+tam_string/np;
				MPI_Send(&msg_ini, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				MPI_Send(&msg_fim, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				pos+=tam_string/np;
				
			}
			else{
				//printf("Ini %d fim %d\n",pos-margem,tam_string-1);
				msg_ini = pos-margem;
				msg_fim = tam_string-1;
				MPI_Send(&msg_ini, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				MPI_Send(&msg_fim, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
				pos+=tam_string/np;
			}
		}
		char buffer[tam_string+1];
		msg_ini = 0;
		msg_fim = (tam_string/np)+margem;
		slice_str(string,buffer,msg_ini,msg_fim);
		result = bmhs(buffer,strlen(buffer),substr,tam_substring);
		printf("\nstring: %s \nmeu pedaco: %s\nsubstr %s rank %d ini %d fim %d result %d\n",string, buffer,substr,meu_rank,msg_ini,msg_fim,result);
		
		fprintf(foutparalelo,"Rank %d\n%s \n%s \n%d \nstring: %s \nmeu pedaco: %s\nsubstr %s rank %d ini %d fim %d result %d\n\n",meu_rank,desc_query,desc_dna,result,string, buffer,substr,meu_rank,msg_ini,msg_fim,result);
		
		fclose(foutparalelo);
	}
	else{
		int msg_ini;
		int msg_fim;
		MPI_Recv(&msg_ini, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&msg_fim, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
		char buffer[tam_string+1];
		slice_str(string,buffer,msg_ini,msg_fim);
		result = bmhs(buffer,strlen(buffer),substr,tam_substring);
		if(result != -1){
			result += msg_ini;
		}
		printf("\nstring: %s \nmeu pedaco: %s\nsubstr %s rank %d ini %d fim %d result %d\n",string, buffer,substr,meu_rank,msg_ini,msg_fim,result);
		fprintf(foutparalelo,"Rank %d\n%s \n%s \n%d \nstring: %s \nmeu pedaco: %s\nsubstr %s rank %d ini %d fim %d result %d\n\n",meu_rank,desc_query,desc_dna,result,string, buffer,substr,meu_rank,msg_ini,msg_fim,result);
		//fprintf(foutparalelo,"%d\n",result);
		fclose(foutparalelo);
	}
	return result;
}

char *bases;
char *str;

int main(int argc, char** argv) {

	bases = (char*) malloc(sizeof(char) * 1000001);
	if (bases == NULL) {
		perror("malloc");
		exit(EXIT_FAILURE);
	}
	str = (char*) malloc(sizeof(char) * 1000001);
	if (str == NULL) {
		perror("malloc str");
		exit(EXIT_FAILURE);
	}

	openfiles();

	char desc_dna[100], desc_query[100];
	char line[100];
	int i, found, result;

	fgets(desc_query, 100, fquery);
	remove_eol(desc_query);
	MPI_Init(&argc, &argv);
	int meu_rank,np;
	MPI_Comm_rank(MPI_COMM_WORLD, &meu_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	char *saida[20];
	sprintf(saida,"dnaparalelo%d.out",meu_rank);
	foutparalelo = fopen(saida, "w");
	fclose(foutparalelo);
	while (!feof(fquery)) {
		fprintf(fout, "%s\n", desc_query);
		// read query string
		fgets(line, 100, fquery);
		remove_eol(line);
		str[0] = 0;
		i = 0;
		do {
			strcat(str + i, line);
			if (fgets(line, 100, fquery) == NULL)
				break;
			remove_eol(line);
			i += 80;
		} while (line[0] != '>');
		strcpy(desc_query, line);
		int meu_rank;
		// read database and search
		found = 0;
		fseek(fdatabase, 0, SEEK_SET);
		fgets(line, 100, fdatabase);
		remove_eol(line);
		while (!feof(fdatabase)) {
			strcpy(desc_dna, line);
			bases[0] = 0;
			i = 0;
			fgets(line, 100, fdatabase);
			remove_eol(line);
			do {
				strcat(bases + i, line);
				if (fgets(line, 100, fdatabase) == NULL)
					break;
				remove_eol(line);
				i += 80;
			} while (line[0] != '>');
			//printf("\nbases\n");
			//printf(bases);
			//printf("\n");
			//printf("str\n");
			//printf(str);
			//result = bmhs(bases, strlen(bases), str, strlen(str));
			result = divide(bases, strlen(bases), str, strlen(str),desc_query,desc_dna);
			if (result > 0) {
				//fflush(stdout);
				//printf("\n%s\nString %s substr %s Rank %d\n%s\n%d\n",desc_query,bases,str,meu_rank,desc_dna, result);
				//fflush(fout);
				fprintf(fout, "%s\n%d\n", desc_dna, result);
				found++;
			}
		}

		if (!found)
			//fflush(stdout);
			//printf("\nRank %d \nNOT FOUND\n",meu_rank);
			//fflush(fout);
			fprintf(fout, "NOT FOUND\n");
	}
	MPI_Finalize();
	
	//if(meu_rank==0){
	//	printf("\nfim mpi\n");
	//	junta_saida(np);
	//}
	closefiles();

	free(str);
	free(bases);

	return EXIT_SUCCESS;
}
