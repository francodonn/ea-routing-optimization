#ifndef INC_REQ_newGA
#define INC_REQ_newGA
#include "newGA.hh"
#include <math.h>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <cstdlib>      // std::rand, std::srand
#include <stdlib.h>     /* srand, rand */

skeleton newGA
{

	// Problem ---------------------------------------------------------------

	Problem::Problem ():_dimension(0)
	{}

	ostream& operator<< (ostream& os, const Problem& pbm)
	{
		os << endl << endl << "Number of Variables " << pbm._dimension
		   << endl;
	
		return os;
	}

	istream& operator>> (istream& is, Problem& pbm)
	{
		char buffer[MAX_BUFFER];
		char buffer2[MAX_BUFFER];
		char buffer3[MAX_BUFFER];
		char buffer4[MAX_BUFFER];
		int i;
		int cantClientes;
		int cantVehiculos;
		char clientsFileName[100];
		char vehiclesInstanceName[100];
		is.getline(buffer,MAX_BUFFER,'\n');
		is.getline(buffer2,MAX_BUFFER,'\n');
		is.getline(buffer3,MAX_BUFFER,'\n');
		is.getline(buffer4,MAX_BUFFER,'\n');
		sscanf(buffer,"%d",&cantClientes);
		sscanf(buffer2,"%d",&cantVehiculos);
		sscanf(buffer3,"%s",clientsFileName);
		sscanf(buffer4,"%s",vehiclesInstanceName);

		pbm._dimension = cantVehiculos * 2 + cantClientes - 1;
		pbm._cant_vehiculos = cantVehiculos;
		pbm._cant_clientes = cantClientes;

		//Pido memoria para almacenar la info de los clientes
		pbm._clients_info = new float *[pbm._cant_clientes];
		for (int i=0;i<pbm._cant_clientes;i++){
			pbm._clients_info[i] = new float [6];
			for (int j=0;j<6;j++)
				pbm._clients_info[i][j]=-1;
		}

		//Pido memoria para almacenar la info de los vehiculos
		pbm._vehicles_info = new float *[pbm._cant_vehiculos];
		for (int i=0; i<pbm._cant_vehiculos; i++){
			pbm._vehicles_info[i] = new float [3];
			for (int j=0;j<3;j++)
				pbm._vehicles_info[i][j]=-1;
		}

		//cargar info vehiculos
		char name [100]= "data/vehiculos/";
		strcat(name, vehiclesInstanceName);
		strcat(name, ".csv");
	    FILE* stream = fopen(name, "r");
		char line[1024];
		float cantVehiculosAgregados = 0;
	    while (fgets(line, 1024, stream))
	    {
			char* tmp = strdup(line);
			const char * vehicle_type = pbm.getfield(tmp,1);
			float vehicle_type_f;
			if(strcmp(vehicle_type, "A")==0) vehicle_type_f = 1;
				else if(strcmp(vehicle_type, "B")==0) vehicle_type_f = 2;
					else if(strcmp(vehicle_type, "C")==0) vehicle_type_f = 3;
						else vehicle_type_f = 0;
			
			free(tmp);
			tmp = strdup(line);
			const char * vehicle_cantidad = pbm.getfield(tmp,2);
			float vehicle_cantidad_f = atof(vehicle_cantidad);

			free(tmp);
			tmp = strdup(line);
			const char * vehicle_capacidad = pbm.getfield(tmp,3);
			float vehicle_capacidad_f = atof(vehicle_capacidad);

			free(tmp);
			tmp = strdup(line);
			const char * vehicle_sueldo = pbm.getfield(tmp,4);
			float vehicle_sueldo_f = atof(vehicle_sueldo);

			for(int i=cantVehiculosAgregados; i<vehicle_cantidad_f + cantVehiculosAgregados; i++){
				pbm._vehicles_info[i][0] = vehicle_type_f;
				pbm._vehicles_info[i][1] = vehicle_capacidad_f;
				pbm._vehicles_info[i][2] = vehicle_sueldo_f;
			}
			cantVehiculosAgregados += vehicle_cantidad_f;

			free(tmp);
		}

		//Cargo el archivo con los datos de clientes
	    char fileName [100]= "data/clients/";
		strcat(fileName, clientsFileName);
		strcat(fileName, ".csv");
	    stream = fopen(fileName, "r");
	    while (fgets(line, 1024, stream))
	    {
	        char* tmp = strdup(line);
	        const char * client_id = pbm.getfield(tmp,1);
	        int client_id_int = atoi(client_id) - 1;

			free(tmp);
			tmp = strdup(line);	
	        const char * coordX= pbm.getfield(tmp,2);
			float coordX_int = (float) atof(coordX);
			pbm._clients_info[client_id_int][0] = coordX_int;

			free(tmp);
			tmp = strdup(line);	
	        const char * coordY= pbm.getfield(tmp,3);
			float coordY_int = atof(coordY);
			pbm._clients_info[client_id_int][1] = coordY_int;

			free(tmp);
			tmp = strdup(line);	
	        const char * demanda= pbm.getfield(tmp,4);
			float demanda_int = atof(demanda);
			pbm._clients_info[client_id_int][2] = demanda_int;

			free(tmp);
			tmp = strdup(line);	
	        const char * inicio= pbm.getfield(tmp,5);
			float inicio_int = atof(inicio);
			pbm._clients_info[client_id_int][3] = inicio_int;

			free(tmp);
			tmp = strdup(line);	
	        const char * fin = pbm.getfield(tmp,6);
			float fin_int = atof(fin);
			pbm._clients_info[client_id_int][4] = fin_int;

			free(tmp);
			tmp = strdup(line);	
	        const char * duracion = pbm.getfield(tmp,7);
			float duracion_int = atof(duracion);
			pbm._clients_info[client_id_int][5] = duracion_int;
			
	        free(tmp);
	    }
		cout<<pbm;


		return is;
	}

	//Funcion para leer el archivo CSV
	const char* Problem::getfield(char* line, int num)
	{
	    const char* tok;
	    for (tok = strtok(line, ",");
	            tok && *tok;
	            tok = strtok(NULL, ",\n"))
	    {
	        if (!--num)
	            return tok;
	    }
	    return NULL;
	}

	bool Problem::operator== (const Problem& pbm) const
	{
		if (_dimension!=pbm.dimension()) return false;
		return true;
	}

	bool Problem::operator!= (const Problem& pbm) const
	{
		return !(*this == pbm);
	}

	Direction Problem::direction() const
	{
		//return maximize;
		return minimize;
	}

	int Problem::dimension() const
	{
		return _dimension;
	}

	int Problem::cant_vehiculos() const
	{
		return _cant_vehiculos;
	}

	int Problem::cant_clientes() const
	{
		return _cant_clientes;
	}

	float ** Problem::clients_info() const
	{
		return _clients_info;
	}

	float ** Problem::vehicles_info() const
	{
		return _vehicles_info;
	}

	Problem::~Problem()
	{
		for(int i=0; i<_cant_clientes; i++){
			delete [] _clients_info[i];
		}
		delete [] _clients_info;
		for(int i=0; i<_cant_vehiculos; i++){
			delete [] _vehicles_info[i];
		}
		delete [] _vehicles_info;
	}

	// Solution --------------------------------------------------------------

	Solution::Solution (const Problem& pbm):_pbm(pbm),_var(pbm.dimension())
	{}

	const Problem& Solution::pbm() const
	{
		return _pbm;
	}

	Solution::Solution(const Solution& sol):_pbm(sol.pbm())
	{
		*this=sol;
	}

	istream& operator>> (istream& is, Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++)
			is >> sol._var[i];
		return is;
	}

	ostream& operator<< (ostream& os, const Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++){
			os << " " << sol._var[i];
		}
		return os;
	}

	NetStream& operator << (NetStream& ns, const Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns << sol._var[i];
		return ns;
	}

	NetStream& operator >> (NetStream& ns, Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns >> sol._var[i];
		return ns;
	}

 	Solution& Solution::operator= (const Solution &sol)
	{
		_var=sol._var;
		return *this;
	}

	bool Solution::operator== (const Solution& sol) const
	{
		if (sol.pbm() != _pbm) return false;
		for(int i = 0; i < _var.size(); i++)
			if(_var[i] != sol._var[i]) return false;
		return true;
	}

	bool Solution::operator!= (const Solution& sol) const
	{
		return !(*this == sol);
	}

	void Solution::initialize()
	{
		
		//Buscamos camion al azar, vemos si entra ese cliente y lo asignamos

		int espacioVehiculo[_pbm.cant_vehiculos()];
		for (int i = 0; i < _pbm.cant_vehiculos(); i++) espacioVehiculo[i] = 0;

		int **  vehiculosAclientes = new int*[_pbm.cant_vehiculos()];
		for (int i = 0; i < _pbm.cant_vehiculos(); i++) {
			vehiculosAclientes[i] = new int[_pbm.cant_clientes()];
			for(int j =0;j<_pbm.cant_clientes();j++) {
				vehiculosAclientes[i][j] = 0;
			}
		}

		for (int i=1; i<_pbm.cant_clientes(); i++){
			int hostVehicle = rand_int(0, _pbm.cant_vehiculos()-1);
			while(espacioVehiculo[hostVehicle]  + _pbm.clients_info()[i][2] > _pbm.vehicles_info()[hostVehicle][1]){
				hostVehicle = rand_int(0, _pbm.cant_vehiculos()-1);
			}
			espacioVehiculo[hostVehicle] += _pbm.clients_info()[i][2];

			vehiculosAclientes[hostVehicle][i] = 1;
		}
		
		int indexVehiculo = 0;
		int i = 0;
		while(indexVehiculo < _pbm.cant_vehiculos()) {
			_var[i] = 0; // el inicio del vehiculo
			i++;
			for(int j=1;j<_pbm.cant_clientes();j++){
				if(vehiculosAclientes[indexVehiculo][j] == 1){
					_var[i] = j;
					i++;
				}
			}
			_var[i] = 0; //el final del vehiculo
			i++;
			indexVehiculo++;
		}

		for (int i = 0; i < _pbm.cant_vehiculos(); i++) {
			delete [] vehiculosAclientes[i];
		}
		delete [] vehiculosAclientes;
	}

	float euclideanDistance(float x1, float x2, float y1, float y2) {
		return sqrt( pow(x1 - x2, 2) + pow(y1 - y2, 2) );
	}

	double Solution::fitness ()
	{
		//Asumimos que las soluciones son todas factibles
        double fitness = 0.0;
		int indexCamionero = 0;
		int i=0;
		bool inicio = true;
		float tiempoTranscurridoVehiculo=0;
		float tiempoEspera=0;
		float penalizacion=0;
		while (i<_var.size()-1){//no consideramos el ultimo nodo ya que es el deposito
			//Ej: [0, 2 ,3,4,0, 0,0 ,0,0 ,0,1,5,0]
			if(_var[i] == 0 && !inicio){
				indexCamionero += 1;
				tiempoTranscurridoVehiculo = 0;
				inicio = true;
			}else{
				if(_var[i] == 0) inicio = false;
				//Multiplicamos distancia con proximo  nodo por sueldo por unidad del conductor
				double distanciaNodos = euclideanDistance(_pbm.clients_info()[_var[i]][0], _pbm.clients_info()[_var[i+1]][0], _pbm.clients_info()[_var[i]][1], _pbm.clients_info()[_var[i+1]][1] ) ;
				tiempoEspera = 0;
				tiempoTranscurridoVehiculo += distanciaNodos;
				if(_pbm.clients_info()[_var[i]][3] > tiempoTranscurridoVehiculo){// si llega antes de que habra el local
					tiempoEspera = _pbm.clients_info()[_var[i]][3] - tiempoTranscurridoVehiculo;
				}
				float tiempoDescarga = _pbm.clients_info()[_var[i]][5];
				tiempoTranscurridoVehiculo += tiempoDescarga + tiempoEspera;

				penalizacion = 0;
				if(_pbm.clients_info()[_var[i]][4] < tiempoTranscurridoVehiculo){ //si llega tarde
					penalizacion = tiempoTranscurridoVehiculo - _pbm.clients_info()[_var[i]][4];
				}

				fitness += (distanciaNodos + tiempoDescarga + tiempoEspera) * _pbm.vehicles_info()[indexCamionero][2] + penalizacion;
			}
			i++;
		}

		return fitness;
	}

	char *Solution::to_String() const
	{
		return (char *)_var.get_first();
	}

	void Solution::to_Solution(char *_string_)
	{
		int *ptr=(int *)_string_;
		for (int i=0;i<_pbm.dimension();i++)
		{
			_var[i]=*ptr;
			ptr++;
		}
	}

	unsigned int Solution::size() const
	{
		return (_pbm.dimension() * sizeof(int));
	}


	int& Solution::var(const int index)
	{
		return _var[index];
	}


	Rarray<int>& Solution::array_var()
	{
		return _var;
	}

	Solution::~Solution()
	{}

	// UserStatistics -------------------------------------------------------

	UserStatistics::UserStatistics ()
	{}

	ostream& operator<< (ostream& os, const UserStatistics& userstat)
	{
		os << "\n---------------------------------------------------------------" << endl;
		os << "                   STATISTICS OF TRIALS                   	 " << endl;
		os << "------------------------------------------------------------------" << endl;

		for (int i=0;i< userstat.result_trials.size();i++)
		{
			os << endl
			   << userstat.result_trials[i].trial
			   << "\t" << userstat.result_trials[i].best_cost_trial
			   << "\t\t" << userstat.result_trials[i].worst_cost_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_evaluation_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_iteration_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].time_best_found_trial
			   << "\t\t" << userstat.result_trials[i].time_spent_trial;
		}
		os << endl << "------------------------------------------------------------------" << endl;
		return os;
	}

	UserStatistics& UserStatistics::operator= (const UserStatistics& userstats)
	{
		result_trials=userstats.result_trials;
		return (*this);
	}

	void UserStatistics::update(const Solver& solver)
	{
		if( (solver.pid()!=0) || (solver.end_trial()!=true)
		  || ((solver.current_iteration()!=solver.setup().nb_evolution_steps())
		       && !terminateQ(solver.pbm(),solver,solver.setup())))
			return;

		struct user_stat *new_stat;

		if ((new_stat=(struct user_stat *)malloc(sizeof(struct user_stat)))==NULL)
			show_message(7);
		new_stat->trial         		 		 = solver.current_trial();
		new_stat->nb_evaluation_best_found_trial = solver.evaluations_best_found_in_trial();
		new_stat->nb_iteration_best_found_trial  = solver.iteration_best_found_in_trial();
		new_stat->worst_cost_trial     		 	 = solver.worst_cost_trial();
		new_stat->best_cost_trial     		 	 = solver.best_cost_trial();
		new_stat->time_best_found_trial		 	 = solver.time_best_found_trial();
		new_stat->time_spent_trial 		 		 = solver.time_spent_trial();

		result_trials.append(*new_stat);
	}

	void UserStatistics::clear()
	{
		result_trials.remove();
	}

	UserStatistics::~UserStatistics()
	{
		result_trials.remove();
	}

// Intra_operator  --------------------------------------------------------------

	Intra_Operator::Intra_Operator(const unsigned int _number_op):_number_operator(_number_op),probability(NULL)
	{}

	unsigned int Intra_Operator::number_operator() const
	{
		return _number_operator;
	}

	Intra_Operator *Intra_Operator::create(const unsigned int _number_op)
	{
		switch (_number_op)
		{
			case 0: return new Crossover;break;
			case 1: return new Mutation();break;
		}
	}

	ostream& operator<< (ostream& os, const Intra_Operator& intra)
	{
		switch (intra.number_operator())
		{
			case 0: os << (Crossover&)intra;break;
			case 1: os << (Mutation&)intra;break;
		}
		return os;
	}

	Intra_Operator::~Intra_Operator()
	{}

//  Crossover:Intra_operator -------------------------------------------------------------

	Crossover::Crossover():Intra_Operator(0)
	{
		probability = new float[1];
	}

	bool check_factibility(Solution sol, Rarray<int>& solArray){
		int indexCamionero = 0;
		int i=0;
		bool recorriendo = false;
	
		int cargaCamioneroActual = 0;

		bool result = true;

		//Obtenemos la carga que lleva cada vehiculo
		while (i < solArray.size()){
			if(solArray[i] == 0 && recorriendo){
				if(cargaCamioneroActual > sol.pbm().vehicles_info()[indexCamionero][1]){
					result = false;
					break;
				}
				indexCamionero += 1;
				recorriendo = false;
				cargaCamioneroActual = 0;
			}else{
				if(solArray[i] == 0) recorriendo  = true;
				cargaCamioneroActual += sol.pbm().clients_info()[solArray[i]][2];
			}
			i++;
		}
		return result;
	}

	void Crossover::cross(Solution& sol1,Solution& sol2) const // dadas dos soluciones de la poblacion, las cruza
	{
		//< "sol2_pre_cross:" << sol2.array_var() << endl;
		int i=0;
		int crossPoint1, crossPoint2;
		int item1, item2, pos1, pos2, posSecond1, posSecond2;
		posSecond2 = 10000;

		//Vamos hasta -2 y empezamos en 1 para evitar el primer y ultimo valor
		crossPoint1 = rand_int(1, sol1.pbm().dimension() - 3);
		crossPoint2 = rand_int(crossPoint1 + 1, sol1.pbm().dimension() - 2);
		
		Rarray<int> aux1( sol1.pbm().dimension()) ;
		Rarray<int> aux2( sol2.pbm().dimension()) ;
		for(int r = 0; r< sol1.pbm().dimension();r++){
			aux1[r] = sol1.var(r); 
			aux2[r] = sol2.var(r); 
		}

		for(int i = crossPoint1; i <= crossPoint2; i++){
			//Get the two items to swap.
			item1 = aux1[i];
			item2 = aux2[i];

			//Intercambiamos haciendo PMX.
			if(item1!=item2 && item1!=0 && item2!=0)
			{
				for(int j = 1; j < aux1.size() - 1; j++){ //menos 2 por el deposito (?)
					if(aux1[j] == item1){
						pos1 = j;
					}else if(aux1[j] == item2){
						pos2 = j;
					}

					if(aux2[j] == item2){
						posSecond1 = j;
					}else if(aux2[j] == item1){
						posSecond2 = j;
					}
				} // j

				aux1[pos1] = item2;

				aux1[pos2] = item1;
			
				aux2[posSecond1] = item1;
				
				aux2[posSecond2] = item2;

			}

		} // i

		if(check_factibility(sol1, aux1) && check_factibility(sol2, aux2) ){
				for(int r = 0; r< sol1.pbm().dimension();r++){
					sol1.var(r) = aux1[r];
					sol2.var(r) = aux2[r];
				}
		}

	}

	void Crossover::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i+1<sols.size();i=i+2)
		 	if (rand01()<=probability[0]) cross(*sols[i],*sols[i+1]);
	}

	ostream& operator<< (ostream& os, const Crossover&  cross)
	{
		 os << "Crossover." << " Probability: "
                    << cross.probability[0]
		    << endl;
		 return os;
	}

	void Crossover::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_crossover_probability",(char *)probability,1,sizeof(float));
	}

	void Crossover::UpdateFromState(const StateCenter& _sc)
	{
		 unsigned long nbytes,length;
		 _sc.get_contents_state_variable("_crossover_probability",(char *)probability,nbytes,length);
	}

	void Crossover::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f ",&op,&probability[0]);
		assert(probability[0]>=0);
	}

	Crossover::~Crossover()
	{
		delete [] probability;
	}

	//  Mutation: Sub_operator -------------------------------------------------------------

	Mutation::Mutation():Intra_Operator(1)
	{
		probability = new float[2];
	}

	void move(Rarray<int>& array, int pos_origin, int pos_dest)
	{
		int item_a_mover = array[pos_origin];
		if(pos_origin < pos_dest) {
			for(int i=pos_origin; i < pos_dest; i++){
				array[i] = array[i+1];
			}
			array[pos_dest] = item_a_mover;
		}
		else{
			for(int i=pos_origin; i > pos_dest; i--){
				array[i] = array[i-1];
			}
			array[pos_dest] = item_a_mover;
		}
	}

	void Mutation::mutate(Solution& sol) const
	{
		//TODO mutacion de ints, check de factibilidad por permutacion
		
		for(int i =0; i<sol.pbm().dimension(); i++){
			if (rand01()<=probability[1]){
				Rarray<int> aux(sol.pbm().dimension());
				for(int r = 0; r< sol.pbm().dimension();r++){
					aux[r] = sol.var(r);
				}
				
				//Hacemos Exchange Mutation
				int index_to_change_origin = rand_int(0, aux.size()-1);
				while(aux[index_to_change_origin]==0){
					index_to_change_origin = rand_int(0, aux.size()-1);
				}

				int index_to_change_dest = rand_int(0, aux.size()-1);
				while(aux[index_to_change_dest]==0 && index_to_change_dest!=index_to_change_origin){
					index_to_change_dest = rand_int(0, aux.size()-1);
				}

				int origin_aux = aux[index_to_change_origin];
				aux[index_to_change_origin]=aux[index_to_change_dest];
				aux[index_to_change_dest] = origin_aux;
				
				//Ahora movemos un cliente de vehiculo a otro
				int index_to_move_origin = rand_int(0, aux.size()-1);
				while(aux[index_to_move_origin]==0){
					index_to_move_origin = rand_int(0, aux.size()-1);
				}

				int index_to_move_dest = rand_int(0, aux.size()-1);
				while(aux[index_to_move_dest]==0 && index_to_move_dest!=index_to_move_origin){
					index_to_move_dest = rand_int(0, aux.size()-1);
				}

				move(sol.array_var(), index_to_move_origin, index_to_move_dest);

				if(check_factibility(sol, aux) ){
					for(int r = 0; r< sol.pbm().dimension();r++){
						sol.var(r) = aux[r];
					}
				}
			}
		}
		
	}

	void Mutation::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i<sols.size();i++)
			if(rand01() <= probability[0])	mutate(*sols[i]);
	}

	ostream& operator<< (ostream& os, const Mutation&  mutation)
	{
		os << "Mutation." << " Probability: " << mutation.probability[0]
		   << " Probability1: " << mutation.probability[1]
		   << endl;
		return os;
	}

	void Mutation::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f %f ",&op,&probability[0],&probability[1]);
		assert(probability[0]>=0);
		assert(probability[1]>=0);
	}

	void Mutation::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_mutation_probability",(char *)probability,2,sizeof(probability));
	}

	void Mutation::UpdateFromState(const StateCenter& _sc)
	{
		unsigned long nbytes,length;
		_sc.get_contents_state_variable("_mutation_probability",(char *)probability,nbytes,length);
	}

	Mutation::~Mutation()
	{
		delete [] probability;
	}

// StopCondition_1 -------------------------------------------------------------------------------------

	StopCondition_1::StopCondition_1():StopCondition()
	{}

	bool StopCondition_1::EvaluateCondition(const Problem& pbm,const Solver& solver,const SetUpParams& setup)
	{
		return false;
	}

	StopCondition_1::~StopCondition_1()
	{}

	//------------------------------------------------------------------------
	// Specific methods ------------------------------------------------------
	//------------------------------------------------------------------------

	bool terminateQ (const Problem& pbm, const Solver& solver,
			 const SetUpParams& setup)
	{
		StopCondition_1 stop;
		return stop.EvaluateCondition(pbm,solver,setup);
	}
}
#endif

