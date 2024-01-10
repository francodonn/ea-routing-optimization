#include "newGA.hh"

int main (int argc, char** argv)
{
	using skeleton newGA;

	system("clear");

	if(argc < 4)
		show_message(1);

	ifstream f1(argv[1]);
	if (!f1) show_message(11);

	ifstream f2(argv[2]);
	if (!f2) show_message(12);

	Problem pbm;
	f2 >> pbm;

	Operator_Pool pool(pbm);
	SetUpParams cfg(pool);
	f1 >> cfg;

	Solver_Seq solver(pbm,cfg);
	solver.run();

	if (solver.pid()==0)
	{
		solver.show_state();
		cout << "Solution: " << solver.global_best_solution()
		     << " Fitness: " << solver.global_best_solution().fitness() << endl;
		cout << "\n\n :( ---------------------- THE END --------------- :) ";

		ofstream fexit(argv[3]);
		if(!fexit) show_message(13);

		Solution sol = solver.global_best_solution();

		//Escribe el resultado en el archivo de salida
		for (int i=0;i<solver.pbm().dimension();i++){
			fexit<< sol.var(i) << " ";
		}
		fexit << endl;
		char vehicleType; 
		int indexCamionero = 0;
		int i=0;
		bool inicio = true;
		int cantClientesVehiculoActual = 0;

		

		while (i<solver.pbm().dimension() ){//no consideramos el ultimo nodo ya que es el deposito
			//Ej: [0, 2 ,3,4,0, 0,0 ,0,0 ,0,1,5,0]
			if(sol.var(i) == 0 && !inicio){
				if(cantClientesVehiculoActual > 0){
					switch((int)solver.pbm().vehicles_info()[indexCamionero][0]){
						case 1: vehicleType = 'A'; break;
						case 2: vehicleType = 'B'; break;
						case 3: vehicleType = 'C'; break;
						default: vehicleType = ' ';
					}
					fexit << vehicleType << " ";
				}
				cantClientesVehiculoActual = 0;
				indexCamionero += 1;
				inicio = true;
			}else{
				if(sol.var(i) == 0) {
					inicio = false;
				}else{
					cantClientesVehiculoActual+=1;
				}
				
				//Multiplicamos distancia con proximo  nodo por sueldo por unidad del conductor
			}
			i++;
		}

		fexit << endl;
		fexit << solver.global_best_solution().fitness() << " ";

	}
	return(0);
}
