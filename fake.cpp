#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <ctime>

int func(int iter);

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));


	for (int i = 0; i < 20; i++)
	{
		std::ofstream terr_in; terr_in.open("terrain.txt");
		int n_hills = 5;//rand()%4 + 1;
		terr_in << n_hills << std::endl;
		for (int p = 0; p < n_hills; p++)
		{
			double r = (double)(rand()%200+50)/1000.0;
			double aprx = (double)(rand()%(100-(int)(2*100*r)) + 100*r);
			double apry = (double)(rand()%(100-(int)(2*100*r)) + 100*r);
			terr_in << aprx/100.0 << " " << apry/100.0 << " " << r << std::endl;
		}
		terr_in.close();

		func(i);
		std::cout << i << std::endl;
	}

	return 0;
}

int func(int iter)
{
	std::ofstream uvalues; uvalues.open("./set/uvalues_" + std::to_string(iter) + ".csv");
	double ps[8] = {0.192, 0.144632, 0.187, 0.2963, 0.0103, 0.021, 0.0164, 0.132368};

	for (int g = 0; g < 8; g++)
	{
		// Parameters.
		int K = 0;
		int Nx_K = 0;
		int Ny_K = 0;
		int hills = 0;
		std::vector<double> xi_x;
		std::vector<double> xi_y;
		std::vector<int>    ii_x;
		std::vector<int>    ii_y;
		std::vector<int>    turb_ids;
		std::vector<double> terrain_dat;
		int Nx, Ny; Nx = Ny = 70;
		double Lx, Ly;
		Lx = Ly = 4.0;
		double dx, dy; dx = Lx / (double)Nx; dy = Ly / (double)Ny;
		double **Utot = new double*[Nx];
		for (int i = 0; i < Nx; i++)
			Utot[i] = new double[Ny]{0.0};
		double **Vtot = new double*[Nx];
		for (int i = 0; i < Nx; i++)
			Vtot[i] = new double[Ny]{0.0};
		double **U0tot = new double*[Nx];
		for (int i = 0; i < Nx; i++)
			U0tot[i] = new double[Ny]{0.0};
		double **V0tot = new double*[Nx];
		for (int i = 0; i < Nx; i++)
			V0tot[i] = new double[Ny]{0.0};
		double **Xtot = new double*[Nx]; 
		for (int i = 0; i < Nx; i++)
			Xtot[i] = new double[Ny]{0.0};
		double **Ytot = new double*[Nx];
		for (int i = 0; i < Nx; i++)
			Ytot[i] = new double[Ny]{0.0};

		// Read xi's.
		std::ifstream turbs; turbs.open("turbines.txt");
		turbs >> Nx_K >> Ny_K >> K;
		// Nx_K = iter+4;
		// Ny_K = iter+4;
		K = Nx_K*Ny_K;
		double **Usij = new double*[K];
		for (int i = 0; i < K; i++)
			Usij[i] = new double[K]{0.0};
		double **Vsij = new double*[K];
		for (int i = 0; i < K; i++)
			Vsij[i] = new double[K]{0.0};

		for (int i = 0; i < K; i++)
		{
			int i_k = i;
			double x_k, y_k;
			x_k = y_k = 0.0;
			// turbs >> i_k;
			turb_ids.push_back(i_k);
			xi_x.push_back(( 0.5+ i_k - (i_k/Nx_K)*Nx_K )*(Lx / (Nx_K+0)));
			xi_y.push_back(( 0.5+ i_k/Nx_K )*(Ly / (Ny_K+0)));
		}

		// Read terrain.
		std::ifstream terrain; terrain.open("terrain.txt");
		terrain >> hills;
		for (int i = 0; i < hills; i++)
		{
			double tk_x, tk_y, size_k;
			tk_x = tk_y = size_k = 0.0;
			terrain >> tk_x >> tk_y >> size_k;
			terrain_dat.push_back(tk_x*Lx);
			terrain_dat.push_back(tk_y*Ly);
			terrain_dat.push_back(size_k);
		}

		// Adjust turbine positions to grid.
		for (int i = 0; i < Nx-1; i++)
		{
			for (int j = 0; j < Ny-1; j++)
			{
				for (int p = 0; p < K; p++)
				{
					if (xi_x[p] >= (i+0.5)*dx && xi_x[p] <= (i+1.5)*dx && xi_y[p] >= (j+0.5)*dy && xi_y[p] <= (j+1.5)*dy)
					{
						xi_x[p] = (i+0.5)*dx;
						ii_x.push_back(i);
						xi_y[p] = (j+0.5)*dy;
						ii_y.push_back(j);
					}
				}
			}
		}

		// Generate field based on terrain.
		std::ofstream out_field; out_field.open("field.txt");
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				double uij, vij;
				uij = vij = 0.0;

				bool isovert = false;
				for (int p = 0; p < hills; p++)
				{
					double xp = (i+0.5)*dx - terrain_dat[3*p + 0];
					double yp = (j+0.5)*dy - terrain_dat[3*p + 1];
					//if (xp*xp + yp*yp > 0.35*pow((terrain_dat[p*hills+2]*sqrt(Lx*Lx + Ly*Ly)),2.0))
					if (xp*xp + yp*yp > pow(terrain_dat[3*p+2]*Lx, 2.0))
					{
						uij += -1.0*(yp) / (xp*xp + yp*yp);
						vij += 1.0*xp / (xp*xp + yp*yp);
					}
					else
						isovert = true;
				}
				Xtot[i][j] = (i+0.5)*dx;
				Ytot[i][j] = (j+0.5)*dy;
				if (!isovert)
				{
					Utot[i][j] = uij + sin(g*M_PI/4.0)*0.75;
					Vtot[i][j] = vij + cos(g*M_PI/4.0)*0.75;
				}
				else
				{
					Utot[i][j] = 0.0;
					Vtot[i][j] = 0.0;
				}
				U0tot[i][j] = Utot[i][j];
				V0tot[i][j] = Vtot[i][j];
			}
		}

		// Generate field based on turbine positions.
		for (int p = 0; p < K; p++)
		{
			for (int i = 0; i < Nx; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					double xp = (i+0.5)*dx - xi_x[p];
					double yp = (j+0.5)*dy - xi_y[p];
					double theta = acos( ((sin(g*M_PI/4.0))*xp + (cos(g*M_PI/4.0))*yp) / ((sqrt(sin(g*M_PI/4.0)*sin(g*M_PI/4.0) + cos(g*M_PI/4.0)*cos(g*M_PI/4.0))) * (sqrt(xp*xp + yp*yp))) );
				
					/*	
					if (xp*xp + yp*yp < 0.75*(sqrt(Lx*Lx+Ly*Ly)) && theta < 40.0*(2*M_PI/360.0) && ((i+0.5)*dx != ii_x[p] && (j+0.5)*dy != ii_y[p]))
					{
						Utot[i][j] = 0.7*(Utot[i][j]);
						Vtot[i][j] = 0.7*(Vtot[i][j]);

						for (int pp = 0; pp < K; pp++)
						{
							if (pp != p && i == ii_x[pp] && j == ii_y[pp])
							{
								std::cout << "C " << p << " " << pp << std::endl;
								Usij[p][pp] = 0.7*U0tot[i][j];
								Vsij[p][pp] = 0.7*V0tot[i][j];
							}
							if (pp == p)
							{
								Usij[p][pp] = U0tot[ii_x[p]][ii_y[p]];
								Vsij[p][pp] = V0tot[ii_x[p]][ii_y[p]];
							}
						}
					}
					*/

					if (xp*xp + yp*yp < 0.75*(sqrt(Lx*Lx + Ly*Ly)) && theta < 45.0*(2.0*M_PI/360.0))
					{
						Utot[i][j] *= 0.7;
						Vtot[i][j] *= 0.7;

						/*
						for (int pp = 0; pp < K; pp++)
						{
							if (p != pp)
							{
								if (xi_x[pp] == xp && xi_y[pp] == yp)
								{
									std::cout << "C: " << p << " " << pp << std::endl;
									//Usij[p][pp] = 0.7*U0tot[i][j];
									//Vsij[p][pp] = 0.7*V0tot[i][j];
								}
							}
							else
							{
								//Usij[p][pp] = 0.0;
								//Vsij[p][pp] = 0.0;
							}
						}
						*/
					}
				}
			}

			for (int pp = 0; pp < K; pp++)
			{
				if (p != pp)
				{
					double xp = xi_x[pp] - xi_x[p];
					double yp = xi_y[pp] - xi_y[p];
					double theta = acos( ((sin(g*M_PI/4.0))*xp + (cos(g*M_PI/4.0))*yp) / ((sqrt(sin(g*M_PI/4.0)*sin(g*M_PI/4.0) + cos(g*M_PI/4.0)*cos(g*M_PI/4.0))) * (sqrt(xp*xp + yp*yp))) );

					if (xp*xp + yp*yp < 0.75*(sqrt(Lx*Lx + Ly*Ly)) && theta < 45.0*(2.0*M_PI/360.0))
					{
						// std::cout << "C: " << p << " " << pp << std::endl;
						Usij[pp][p] = 0.7*U0tot[ii_x[pp]][ii_y[pp]];
						Vsij[pp][p] = 0.7*V0tot[ii_x[pp]][ii_y[pp]];
					}
				}
				else
				{
					Usij[p][p] = U0tot[ii_x[p]][ii_y[p]];
					Vsij[p][p] = V0tot[ii_x[p]][ii_y[p]];
				}
			}

			// std::cout << xi_x[p] << " " << xi_y[p] << std::endl;
		}

		// Print.
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				out_field << (i+0.5)*dx << " " << (j+0.5)*dy << " " << (U0tot[i][j] / 5.0)*dx << " " << (V0tot[i][j] / 5.0)*dy << " " << sqrt(U0tot[i][j]*U0tot[i][j] + V0tot[i][j]*V0tot[i][j]) << std::endl;
			}
			out_field << std::endl;
		}
		out_field.close();
		std::cout << "Check..." << std::endl;
		system("gnuplot.exe plot.p");
		std::string sysplot = "cp plot.png ./set/plot_" + std::to_string(iter) + ".png";
		std::string syscopy = "cp terrain.txt ./set/terrain_" + std::to_string(iter) + ".txt";
		system(sysplot.c_str());
		system(syscopy.c_str());

		// uvalues << ps[g] << "\n";
		for (int p = 0; p < K; p++)
		{
			int ip = ii_x[p];
			int jp = ii_y[p];
			for (int pp = 0; pp < K; pp++)
			{
				if (Usij[p][pp] == 0.0 && Vsij[p][pp] == 0.0)
				{
					Usij[p][pp] = U0tot[ip][jp];
					Vsij[p][pp] = V0tot[ip][jp];
				}
				//uvalues << sqrt(U0tot[ip][jp]*U0tot[ip][jp] + V0tot[ip][jp]*V0tot[ip][jp]) << "," << sqrt(Usij[p][pp]*Usij[p][pp] + Vsij[p][pp]*Vsij[p][pp]) << " ";
			}
			//uvalues << std::endl;
		}
		//uvalues << std::endl;
		
		for (int pp = 0; pp < Nx_K*Ny_K; pp++)
		{
			for (int p = 0; p < Nx_K*Ny_K; p++)
			{
				uvalues << ps[g] << ",";
				std::vector<int>::iterator it_pp = std::find(turb_ids.begin(), turb_ids.end(), pp);
				std::vector<int>::iterator it_p = std::find(turb_ids.begin(), turb_ids.end(), p);
				if (it_pp != turb_ids.end())
				{
					int id_pp = it_pp - turb_ids.begin();
					int ip = ii_x[id_pp];
					int jp = ii_y[id_pp];

					if (it_p != turb_ids.end())
					{
						int id_p = it_p - turb_ids.begin();

						uvalues << sqrt(U0tot[ip][jp]*U0tot[ip][jp] + V0tot[ip][jp]*V0tot[ip][jp]) << "," << sqrt(Usij[id_pp][id_p]*Usij[id_pp][id_p] + Vsij[id_pp][id_p]*Vsij[id_pp][id_p]) << std::endl;

					}
					else
					{
						uvalues << sqrt(U0tot[ip][jp]*U0tot[ip][jp] + V0tot[ip][jp]*V0tot[ip][jp]) << ",0" << std::endl;
					}
				}
				else
				{
					uvalues << "0,0" << std::endl;
				}
			}
		}
	}
	uvalues.close();

	return 0;
}
