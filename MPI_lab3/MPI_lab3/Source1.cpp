#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <ctime>
using namespace std;

#define _CRT_SECURE_NO_WARNINGS

void Kern(double kernel[3][3], int radius, double sigma)
{
	double norm = 0;

	for (int i = -radius; i <= radius; i++)
		for (int j = -radius; j <= radius; j++)
		{
			kernel[i + radius][j + radius] = (exp(-(i * i + j * j) / (2 * sigma * sigma)));
			norm += kernel[i + radius][j + radius];
		}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			kernel[i][j] /= norm;
}

int Clamp(int value, int min, int max)
{
	if (value < min)
		return min;

	if (value > max)
		return max;

	return value;
}

int main(int argc, char* argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double sigma = 1.0;
	int diam = 3;
	int q = 20;
	int rows = q, cols = q;
	int *image = (int*)malloc(rows*cols * sizeof(int*));
	int *result2 = (int*)malloc(rows*cols * sizeof(int*));
	int **img;
	double kernel[3][3];


	if (rank == 0)
	{
		img = (int**)malloc(rows * sizeof(int*));
		for (int i = 0; i < rows; ++i)
			img[i] = (int*)malloc(cols * sizeof(int));
		srand(time(NULL));
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				img[i][j] = rand() % 256;
			}

		int k = 0;

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				image[k] = img[i][j];
				k++;
			}

		Kern(kernel, 1, sigma);

	}

	if (rank == 0)
	{
		double temp;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				temp = 0;
				for (int q = -1; q <= 1; q++)
					for (int l = -1; l <= 1; l++)
					{
						int idX = Clamp(i + q, 0, rows - 1);
						int idY = Clamp(j + l, 0, cols - 1);
						temp += image[idX*rows + idY] * kernel[q + 1][l + 1];
					}
				result2[i*rows + j] = int(temp);
			}
	}

	MPI_Bcast(kernel, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Status status;
	MPI_Datatype type, type2;

	int count = ((cols / size) + 2);

	MPI_Type_contiguous(count, MPI_INT, &type2);
	MPI_Type_commit(&type2);

	int stride = sizeof(int) * cols;

	MPI_Type_create_hvector(rows, 1, stride, type2, &type);
	MPI_Type_commit(&type);

	int tmp = cols - ((cols / size)*(size - 2) + (cols / size + 1)) - 1;

	double start = MPI_Wtime();

	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Send((image + ((tmp - i + 1) + (i - 1)*(cols / size + 1))), 1, type, i, NULL, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(image, 1, type, 0, NULL, MPI_COMM_WORLD, &status);
	}

	/*for (int i = 1; i < size; i++)
	{
	MPI_Barrier(MPI_COMM_WORLD);
	if (i == rank)
	{
	cout << i << endl;
	for (int k = 0; k < rows*(cols); k++)
	{
	cout << image[k] << " ";
	}
	cout << endl;
	}
	}*/
	int temp = count;
	if (rank == 0)
	{
		temp = tmp+2;
	}

	int *result = (int*)malloc(rows*temp * sizeof(int*));
	/*for (int q = 0; q < size; q++)
	{
	MPI_Barrier(MPI_COMM_WORLD);
	if (q == rank)
	{
	cout << q << endl;
	for (int k = -1; k <= 1; k++)
	for (int q = -1; q <= 1; q++)
	{
	cout << kernel[1 + k][1 + q] << endl;
	}
	//cout << q << endl;
	}
	}
	*/

	double qwerty;
	//if (rank == 3)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < temp; j++)
			{
				//cout << "            " << i*rows + j << endl;
				qwerty = 0;
				for (int k = -1; k <= 1; k++)
					for (int q = -1; q <= 1; q++)
					{
						int idX = Clamp(i + k, 0, rows - 1);
						int idY = Clamp(j + q, 0, temp - 1);

						qwerty += image[idX*rows + idY] * kernel[1 + k][1 + q];

						//cout << idX << " " << idY << endl;
						//cout << "      " <<image[idX*rows + idY] << endl;
					}
				result[i*temp + j] = Clamp(qwerty, 0, 255);
				//cout << result[i*temp + j] << endl;
			}
	}


	/*MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < temp; j++)
			{
				std::cout << result[i*temp + j] << "  ";
			}
			std::cout << endl;
		}
	}*/
	//MPI_Barrier(MPI_COMM_WORLD);


	if (rank != 0)
	{
		MPI_Send(result, rows*count, MPI_INT, 0, NULL, MPI_COMM_WORLD);
	}
	if (rank == 0)
	{

		/*cout << temp;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < temp; j++)
			{
				image[i*rows + j] = result[i*temp + j];
			}*/
		temp--;
		int *res = (int*)malloc(rows*count * sizeof(int*));
		for (int s = 1; s < size; s++)
		{

			MPI_Recv(res, rows*count, MPI_INT, s, NULL, MPI_COMM_WORLD, &status);

			/*if (s == 1)
			{
				for (int i = 0; i < rows; i++)
				{
					image[i*rows + ((j)+temp + (count - 2)*((s - 1)))] = result[i*count + j + 1];
				}
			}*/
			if (s == size-1)
			{
				//count--;
			}
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < count; j++)
				{
					image[i*rows + ((j)+temp + (count - 2)*((s - 1)))] = res[i*count + j + 1];
				}
			}
					for (int i = 0; i < rows; i++)
			for (int j = 0; j < temp; j++)
			{
				image[i*rows + j] = result[i*temp + j];
			}

			//temp -= 2;
		}
		temp++;
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < temp-1; j++)
			{
				image[i*rows + j] = result[i*temp + j];
			}

		printf("time = %f\n", (MPI_Wtime() - start));
	}

	if (rank == 0)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				if (result2[i*rows + j] - image[i*rows + j] != 0)
				{
					cout << "       " << i << "  " << j << endl;
					//cout << result2[i*rows + j] << endl;
					//cout << image[i*rows + j] << endl;

					cout << "false" << endl;
				}
			}
	}


	/*if (rank == 0)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				std::cout << image[i* rows + j] << "  ";
			}
			std::cout << endl;
		}

		std::cout << endl;
		std::cout << endl;
		std::cout << endl;
	}


	if (rank == 0)
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				std::cout << result2[i* rows + j] << "  ";
			}
			std::cout << endl;
		}
	}*/
	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}