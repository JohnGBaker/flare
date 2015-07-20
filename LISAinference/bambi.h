#ifndef __BAMBI_H__
#define __BAMBI_H__ 1

void quickSort(float numbers[], int array_size);
void q_sort(float numbers[], int left, int right);
void FirstRunCheck(int);
void GetNetTol();
//void CopyFile(std::string, std::string);
void CopyFile(char *, char *);

void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context);
void bambi(int *ndata, int *ndim, double **BAMBIData, double *lowlike);

#endif
