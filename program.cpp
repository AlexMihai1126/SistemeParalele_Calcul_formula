#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>

#define N 200
#define TAG_VECTOR_X 1
#define TAG_VECTOR_Y 2
#define TAG_MATRIX_A 3
#define MASTER_PROC_ID 0

using namespace std;
// serializam matricea rowsToSend linii si N coloane
// Serialize a portion of the matrix into a 1D array, starting from a given offset row
void serialize(const vector<vector<float>> &data, float *buffer, int offset, int rowsToSend)
{
    int index = 0;
    for (int i = offset; i < offset + rowsToSend; ++i)
    {
        for (float value : data[i])
        {
            buffer[index++] = value;
        }
    }
}

// Deserialize a 1D array into a vector<vector<float>> matrix
void deserialize(const float *buffer, vector<vector<float>> &data, int rowsToReceive)
{
    data.clear();
    int index = 0;
    for (int i = 0; i < rowsToReceive; ++i)
    {
        vector<float> row;
        for (int j = 0; j < N; ++j)
        {
            row.push_back(buffer[index++]);
        }
        data.push_back(row);
    }
}

int main(int argc, char **argv)
{
    int mainRank, nrProc, secondRank;
    MPI_Comm lowComm, highComm;

    // init mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mainRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrProc);

    if (nrProc % 2 != 0 || nrProc < 4)
    {
        cout << "Numarul de procese trebuie sa fie par si mai mare ca 2!" << endl;
        MPI_Finalize();
        return 0;
    } // verificam daca nr de procese este par inainte de a face alocari
    // verificam si ca numarul minim de procese sa fie 4 deoarece neavand voie ca procesul master pe fiecare subgrup
    // sa faca calcule ci doar sa adune si sa trimita date mai departe, am ajunge la o eroare unde nu avem procese worker.

    int lowGroupMaxRank = nrProc / 2;
    int highGroupMinRank = nrProc - lowGroupMaxRank;

    if (mainRank < lowGroupMaxRank)
    {
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &lowComm);
        MPI_Comm_rank(lowComm, &secondRank);
    }
    else
    {
        MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &highComm);
        MPI_Comm_rank(highComm, &secondRank);
    }

    if (mainRank < lowGroupMaxRank)
    {
        // cout<<"Low group process with rank "<<mainRank<<" and low group rank "<<secondRank<<endl;
        int masterLow = 0;
        int elementsPerProcess = N / (lowGroupMaxRank - 1); // calculam numarul de elemente care va fi procesat de fiecare worker
        int remainder = N % (lowGroupMaxRank - 1);
        if (mainRank == lowGroupMaxRank - 1)
        {
            elementsPerProcess += remainder;
        } // daca nu se poate imparti egal, ultimul proces va primi restul de elemente

        if (secondRank == 0 && mainRank < lowGroupMaxRank)
        {
            int masterLowGroup = mainRank;
            vector<float> xData;
            vector<float> yData;

            ifstream xFile("x.dat");
            if (!xFile.is_open())
            {
                cerr << "Error: Unable to open x.dat for reading." << endl;
                MPI_Finalize();
                return 0;
            }
            float value;
            while (xFile >> value)
            {
                xData.push_back(value);
            }
            xFile.close();

            ifstream yFile("y.dat");
            if (!yFile.is_open())
            {
                cerr << "Error: Unable to open y.dat for reading." << endl;
                MPI_Finalize();
                return 0;
            }
            while (yFile >> value)
            {
                yData.push_back(value);
            }
            yFile.close();

            // Send data to every process in the lowComm group
            for (int destRank = 1; destRank < lowGroupMaxRank; ++destRank) // Start from 1 to skip sending to itself
            {
                // Calculate the offset for the current destination process
                int offset = (destRank - 1) * elementsPerProcess;

                // Send xData
                MPI_Send(xData.data() + offset, elementsPerProcess, MPI_FLOAT, destRank, TAG_VECTOR_X, lowComm);

                // Send yData
                MPI_Send(yData.data() + offset, elementsPerProcess, MPI_FLOAT, destRank, TAG_VECTOR_Y, lowComm);
            }
        }
        else
        {
            // Allocate memory for received data
            vector<float> localXData(elementsPerProcess);
            vector<float> localYData(elementsPerProcess);

            // Receive xData
            MPI_Recv(localXData.data(), elementsPerProcess, MPI_FLOAT, masterLow, TAG_VECTOR_X, lowComm, MPI_STATUS_IGNORE);

            // Receive yData
            MPI_Recv(localYData.data(), elementsPerProcess, MPI_FLOAT, masterLow, TAG_VECTOR_Y, lowComm, MPI_STATUS_IGNORE);

            /* // Process received data
            cout << "Process " << mainRank << " received " << elementsPerProcess << " elements of xData and yData."
                 << endl;
            for (int i = 0; i < elementsPerProcess; ++i)
            {
                cout << "xData[" << i << "]: " << localXData[i] << ", yData[" << i << "]: " << localYData[i] << endl;
            } */
        }
    }
    else
    {
        if (mainRank >= highGroupMinRank)
        {
            // cout<<"High group process with rank "<<mainRank<<" and high group rank "<<secondRank<<endl;
            int masterHigh = 0;
            int elementsPerProcess = N / ((nrProc - highGroupMinRank) - 1); // calculam numarul de elemente care va fi procesat de fiecare worker
            int remainder = N % ((nrProc - highGroupMinRank) - 1);
            if (mainRank == (nrProc - 1))
            {
                elementsPerProcess += remainder;
            } // daca nu se poate imparti egal, ultimul proces va primi restul de elemente

            if (secondRank == 0 && mainRank >= highGroupMinRank)
            {
                vector<float> xData;
                vector<float> yData;
                vector<vector<float>> matData(N, vector<float>(N));

                ifstream xFile("x.dat");
                if (!xFile.is_open())
                {
                    cerr << "Error: Unable to open x.dat for reading." << endl;
                    MPI_Finalize();
                    return 0;
                }
                float value;
                while (xFile >> value)
                {
                    xData.push_back(value);
                }
                xFile.close();

                ifstream yFile("y.dat");
                if (!yFile.is_open())
                {
                    cerr << "Error: Unable to open y.dat for reading." << endl;
                    MPI_Finalize();
                    return 0;
                }
                while (yFile >> value)
                {
                    yData.push_back(value);
                }
                yFile.close();

                ifstream file("mat.dat");
                if (!file.is_open())
                {
                    cerr << "Error: Unable to open mat.dat file." << endl;
                    MPI_Finalize();
                    return 0;
                }

                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        file >> matData[i][j];
                    }
                }
                file.close();

                for (int destRank = 1; destRank < lowGroupMaxRank; ++destRank)
                {
                    int offset = (destRank - 1) * elementsPerProcess;
                    MPI_Send(xData.data() + offset, elementsPerProcess, MPI_FLOAT, destRank, TAG_VECTOR_X, highComm);
                    MPI_Send(yData.data(), N, MPI_FLOAT, destRank, TAG_VECTOR_Y, highComm);

                    // Serialize matrix data, starting from offset row
                    auto *buffer = new float[elementsPerProcess * N];
                    serialize(matData, buffer, offset, elementsPerProcess);
                    MPI_Send(buffer, elementsPerProcess * N, MPI_FLOAT, destRank, TAG_MATRIX_A, highComm);
                    delete[] buffer;
                }
            }
            else
            {
                // primim datele
                vector<float> localXData(elementsPerProcess);
                vector<float> localYData(N);
                vector<vector<float>> localMatData(elementsPerProcess, vector<float>(N));
                
                MPI_Recv(localXData.data(), elementsPerProcess, MPI_FLOAT, masterHigh, TAG_VECTOR_X, highComm, MPI_STATUS_IGNORE);
                MPI_Recv(localYData.data(), N, MPI_FLOAT, masterHigh, TAG_VECTOR_Y, highComm, MPI_STATUS_IGNORE);

                // Allocate memory for the buffer
                auto *buffer = new float[elementsPerProcess * N];
                MPI_Recv(buffer, elementsPerProcess * N, MPI_FLOAT, masterHigh, TAG_MATRIX_A, highComm, MPI_STATUS_IGNORE);

                // Deserialize matrix data
                deserialize(buffer, localMatData, elementsPerProcess);
                delete[] buffer;

                // Print received data
                /* cout << "Process " << mainRank << " received X vector:" << endl;
                for (int i = 0; i < elementsPerProcess; ++i)
                {
                    cout << localXData[i] << " ";
                }
                cout << endl;

                cout << "Process " << mainRank << " received Y vector:" << endl;
                for (int i = 0; i < N; ++i)
                {
                    cout << localYData[i] << " ";
                }
                cout << endl; */

                /* cout << "Process " << mainRank << " received Matrix data:" << endl;
                for (int i = 0; i < elementsPerProcess; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        cout << localMatData[i][j] << " ";
                    }
                    cout << endl;
                } */
            
        }
    }
    MPI_Finalize();

    return 0;
}
