/*
 * Simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.2
 *
 * (c) 2020 Arturo Gonzalez Escribano
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <cputils.h>
#include <mpi.h>

/* function use for time data output */
#if !defined(CP_TABLON)
//define ansi color for better output
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"
void printTimes(double *timeGatherer, int nprocs, double ttotal)
{
	double totalGatheredTime = 0;
	int i = 0;
	for (i = 0; i < nprocs; i++)
	{
		printf("%0.6lf\t", timeGatherer[i]);
		totalGatheredTime += timeGatherer[i];
	}
	printf("\n\t");
	for (i = 0; i < nprocs; i++)
		printf("%0.6lf%%\t", timeGatherer[i] / ttotal * 100);
	printf("\n\tTotal: \t\t" ANSI_COLOR_GREEN "%0.6lf" ANSI_COLOR_RESET " \n\tPorcentaje:\t" ANSI_COLOR_GREEN "%0.6lf%%" ANSI_COLOR_RESET "\n", totalGatheredTime / nprocs, (totalGatheredTime / nprocs) / ttotal * 100);
}
#endif
/* Structure to store data of a cell */
typedef struct
{
	float pos_row, pos_col;		  // Position
	float mov_row, mov_col;		  // Direction of movement
	float choose_mov[3];		  // Genes: Probabilities of 0 turning-left; 1 advance; 2 turning-right
	float storage;				  // Food/Energy stored
	int age;					  // Number of steps that the cell has been alive
	unsigned short random_seq[3]; // Status value of its particular random sequence
	bool alive;					  // Flag indicating if the cell is still alive
} Cell;
/* prototype function to compare cells by position*/
int compareCells(const void *_a, const void *_b);

/* Structure for simulation statistics */
typedef struct
{
	int history_total_cells;	 // Accumulated number of cells created
	int history_dead_cells;		 // Accumulated number of dead cells
	int history_max_alive_cells; // Maximum number of cells alive in a step
	int history_max_new_cells;	 // Maximum number of cells created in a step
	int history_max_dead_cells;	 // Maximum number of cells died in a step
	int history_max_age;		 // Maximum age achieved by a cell
	float history_max_food;		 // Maximum food level in a position of the culture
} Statistics;

/*
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be changed and/or optimized by the students
 *
 */
#define accessMat(arr, exp1, exp2) arr[(int)(exp1)*columns + (int)(exp2)]

/*
 * Function: Choose a new direction of movement for a cell
 * 	This function can be changed and/or optimized by the students
 */
void cell_new_direction(Cell *cell)
{
	float angle = (float)(2 * M_PI * erand48(cell->random_seq));
	cell->mov_row = sinf(angle);
	cell->mov_col = cosf(angle);
}

/*
 * Function: Mutation of the movement genes on a new cell
 * 	This function can be changed and/or optimized by the students
 */
void cell_mutation(Cell *cell)
{
	/* 1. Select which genes change:
	 	0 Left grows taking part of the Advance part
	 	1 Advance grows taking part of the Left part
	 	2 Advance grows taking part of the Right part
	 	3 Right grows taking part of the Advance part
	*/
	int mutation_type = (int)(4 * erand48(cell->random_seq));
	/* 2. Select the amount of mutation (up to 50%) */
	float mutation_percentage = (float)(0.5 * erand48(cell->random_seq));
	/* 3. Apply the mutation */
	float mutation_value;
	switch (mutation_type)
	{
	case 0:
		mutation_value = cell->choose_mov[1] * mutation_percentage;
		cell->choose_mov[1] -= mutation_value;
		cell->choose_mov[0] += mutation_value;
		break;
	case 1:
		mutation_value = cell->choose_mov[0] * mutation_percentage;
		cell->choose_mov[0] -= mutation_value;
		cell->choose_mov[1] += mutation_value;
		break;
	case 2:
		mutation_value = cell->choose_mov[2] * mutation_percentage;
		cell->choose_mov[2] -= mutation_value;
		cell->choose_mov[1] += mutation_value;
		break;
	case 3:
		mutation_value = cell->choose_mov[1] * mutation_percentage;
		cell->choose_mov[1] -= mutation_value;
		cell->choose_mov[2] += mutation_value;
		break;
	default:
		fprintf(stderr, "Error: Imposible type of mutation\n");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	/* 4. Correct potential precision problems */
	cell->choose_mov[2] = 1.0f - cell->choose_mov[1] - cell->choose_mov[0];
}

#ifdef DEBUG
/*
 * Function: Print the current state of the simulation
 */
void print_status(int iteration, int rows, int columns, float *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat)
{
	/*
	 * You don't need to optimize this function, it is only for pretty printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i, j;

	printf("Iteration: %d\n", iteration);
	printf("+");
	for (j = 0; j < columns; j++)
		printf("---");
	printf("+\n");
	for (i = 0; i < rows; i++)
	{
		printf("|");
		for (j = 0; j < columns; j++)
		{
			char symbol;
			if (accessMat(culture, i, j) >= 20)
				symbol = '+';
			else if (accessMat(culture, i, j) >= 10)
				symbol = '*';
			else if (accessMat(culture, i, j) >= 5)
				symbol = '.';
			else
				symbol = ' ';

			int t;
			int counter = 0;
			for (t = 0; t < num_cells; t++)
			{
				int row = (int)(cells[t].pos_row);
				int col = (int)(cells[t].pos_col);
				if (cells[t].alive && row == i && col == j)
				{
					counter++;
				}
			}
			if (counter > 9)
				printf("(M)");
			else if (counter > 0)
				printf("(%1d)", counter);
			else
				printf(" %c ", symbol);
		}
		printf("|\n");
	}
	printf("+");
	for (j = 0; j < columns; j++)
		printf("---");
	printf("+\n");
	printf("Num_cells_alive: %04d\nHistory( Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f )\n\n",
		   num_cells_alive,
		   sim_stat.history_total_cells,
		   sim_stat.history_dead_cells,
		   sim_stat.history_max_alive_cells,
		   sim_stat.history_max_new_cells,
		   sim_stat.history_max_dead_cells,
		   sim_stat.history_max_age,
		   sim_stat.history_max_food);
}
#endif

/*
 * Function: Print usage line in stderr
 */
void show_usage(char *program_name)
{
	fprintf(stderr, "Usage: %s ", program_name);
	fprintf(stderr, "<rows> <columns> <maxIter> <max_food> <food_density> <food_level> <short_rnd1> <short_rnd2> <short_rnd3> <num_cells>\n");
	fprintf(stderr, "\tOptional arguments for special food spot: [ <row> <col> <size_rows> <size_cols> <density> <level> ]\n");
	fprintf(stderr, "\n");
}

/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[])
{
	int i, j;

	// Simulation data
	int max_iter;		  // Maximum number of simulation steps
	int rows, columns;	  // Cultivation area sizes
	float *culture;		  // Cultivation area values
	short *culture_cells; // Ancillary structure to count the number of cells in a culture space

	float max_food;		// Maximum level of food on any position
	float food_density; // Number of food sources introduced per step
	float food_level;	// Maximum number of food level in a new source

	bool food_spot_active = false;	// Special food spot: Active
	int food_spot_row = 0;			// Special food spot: Initial row
	int food_spot_col = 0;			// Special food spot: Initial row
	int food_spot_size_rows = 0;	// Special food spot: Rows size
	int food_spot_size_cols = 0;	// Special food spot: Cols size
	float food_spot_density = 0.0f; // Special food spot: Food density
	float food_spot_level = 0.0f;	// Special food spot: Food level

	unsigned short init_random_seq[3];		// Status of the init random sequence
	unsigned short food_random_seq[3];		// Status of the food random sequence
	unsigned short food_spot_random_seq[3]; // Status of the special food spot random sequence

	int num_cells; // Number of cells currently stored in the list
	Cell *cells;   // List to store cells information

	// Data needed for parallel task
	int tag;
	float *my_culture;			// Cultivation area values PER PROCESS
	short *my_culture_cells;		// Ancillary structure to count the number of cells in a culture space PER PROCESS
	int my_rows_number = 0;
	int rows_left = 0;
	int my_begin = 0;
	int my_end;

	// Statistics
	Statistics sim_stat;
	sim_stat.history_total_cells = 0;
	sim_stat.history_dead_cells = 0;
	sim_stat.history_max_alive_cells = 0;
	sim_stat.history_max_new_cells = 0;
	sim_stat.history_max_dead_cells = 0;
	sim_stat.history_max_age = 0;
	sim_stat.history_max_food = 0.0f;

	/* 0. Initialize MPI */
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* 1. Read simulation arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 11)
	{
		fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage(argv[0]);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	/* 1.2. Read culture sizes, maximum number of iterations */
	rows = atoi(argv[1]);
	columns = atoi(argv[2]);
	max_iter = atoi(argv[3]);

	/* 1.3. Food data */
	max_food = atof(argv[4]);
	food_density = atof(argv[5]);
	food_level = atof(argv[6]);

	/* 1.4. Read random sequences initializer */
	for (i = 0; i < 3; i++)
	{
		init_random_seq[i] = (unsigned short)atoi(argv[7 + i]);
	}

	/* 1.5. Read number of cells */
	num_cells = atoi(argv[10]);

	/* 1.6. Read special food spot */
	if (argc > 11)
	{
		if (argc < 17)
		{
			fprintf(stderr, "-- Error in number of special-food-spot arguments in the command line\n\n");
			show_usage(argv[0]);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		else
		{
			food_spot_active = true;
			food_spot_row = atoi(argv[11]);
			food_spot_col = atoi(argv[12]);
			food_spot_size_rows = atoi(argv[13]);
			food_spot_size_cols = atoi(argv[14]);
			food_spot_density = atof(argv[15]);
			food_spot_level = atof(argv[16]);

			// Check non-used trailing arguments
			if (argc > 17)
			{
				fprintf(stderr, "-- Error: too many arguments in the command line\n\n");
				show_usage(argv[0]);
				MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
			}
		}
	}

#ifdef DEBUG
	/* 1.7. Print arguments */
	printf("Arguments, Rows: %d, Columns: %d, max_iter: %d\n", rows, columns, max_iter);
	printf("Arguments, Max.food: %f, Food density: %f, Food level: %f\n", max_food, food_density, food_level);
	printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", init_random_seq[0], init_random_seq[1], init_random_seq[2]);
	if (food_spot_active)
	{
		printf("Arguments, Food_spot, pos(%d,%d), size(%d,%d), Density: %f, Level: %f\n",
			   food_spot_row, food_spot_col, food_spot_size_rows, food_spot_size_cols, food_spot_density, food_spot_level);
	}
	printf("Initial cells: %d\n", num_cells);
#endif // DEBUG

	/* 1.8. Initialize random sequences for food dropping */
	for (i = 0; i < 3; i++)
	{
		food_random_seq[i] = (unsigned short)nrand48(init_random_seq);
		food_spot_random_seq[i] = (unsigned short)nrand48(init_random_seq);
	}

	/* 1.9. Initialize random sequences of cells */
	cells = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
	if (cells == NULL)
	{
		fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
	for (i = 0; i < num_cells; i++)
	{
		// Initialize the cell ramdom sequences
		for (j = 0; j < 3; j++)
			cells[i].random_seq[j] = (unsigned short)nrand48(init_random_seq);
	}

#ifdef DEBUG
	/* 1.10. Print random seed of the initial cells */
	/*
	printf("Initial cells random seeds: %d\n", num_cells );
	for( i=0; i<num_cells; i++ )
		printf("\tCell %d, Random seq: %hu,%hu,%hu\n", i, cells[i].random_seq[0], cells[i].random_seq[1], cells[i].random_seq[2] );
	*/
#endif // DEBUG

	/* 2. Start global timer */
	MPI_Barrier(MPI_COMM_WORLD);
	double ttotal = cp_Wtime();

	/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */

	//TODO: Balance the number of cells in each processor
	//TODO: Split culture and culture cells
	//TODO: Make proccessor 0 the manager of food number generation and statistics collection (corner case: only one proccessor)

#if !defined(CP_TABLON)	  //Precompilaciñon para evitar medir tiempos en el momento de la prueba en el servidor http://frontendv.infor.uva.es/faq#6 \
						  // 2.1 Time variables for single loop iterations
	double timeInitCS;	  // 3.1 Initialize culture surface
	double timeInitCells; // 3.2 Initialize cells

	double timeML;				  // 1º Loop of the Simultation. Main Loop
	double timeNormalSpreadingL;  // NormalSpreading Loop
	double timeSpecialSpreadingL; // SpecialSpreading Loop
	double timeCellMovementL;	  // CellMovement Loop
	double timeCellActionsL;	  // CellActions Loop
	double timeMovingAliveCellsL; // MovingAliveCells Loop
	double timeJoinCellsListL;	  // JoinCellsList Loop
	double timeDecreaseFoodL;	  // DecreaseFood Loop
	double timeDataCollectionL;	  // Sync statistics

	// 2.2 Time variables for total time invested in each loop

	//THESE DON'T NEED A SPECIAL COUNTING
	//timeInitCS;				// 3.1 Initialize culture surface
	//timeInitCells;		// 3.2 Initialize cells
	//timeML = 0.0f;  // 1º Loop of the Simultation. Main Loop

	double timeNormalSpreadingT = 0.0f;	 // NormalSpreading Loop
	double timeSpecialSpreadingT = 0.0f; // SpecialSpreading Loop
	double timeCellMovementT = 0.0f;	 // CellMovement Loop
	double timeCellActionsT = 0.0f;		 // CellActions Loop
	double timeMovingAliveCellsT = 0.0f; // MovingAliveCells Loop
	double timeJoinCellsListT = 0.0f;	 // JoinCellsList Loop
	double timeDecreaseFoodT = 0.0f;	 // DecreaseFood Loop
	double timeDataCollectionT = 0.0f;	 // Sync statistics

#endif

	//Variable for general check of cell alive
	bool any_cell_alive = false;

	if (num_cells > 0)
	{
		any_cell_alive = true;
	}
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	/* Create MPI_Cell datatype */
	// Number of field blocks
	int fields = 4;
	// Number of elements per block
	int array_of_blocklengths[] = {8, 1, 3, 1};
	// Block displacements
	MPI_Aint array_of_displacements[] = {
		offsetof(Cell, pos_row),
		offsetof(Cell, age),
		offsetof(Cell, random_seq),
		offsetof(Cell, alive),
	};
	// Block types
	MPI_Datatype array_of_types[] = {MPI_FLOAT, MPI_INT, MPI_UNSIGNED_SHORT, MPI_C_BOOL};
	MPI_Datatype MPI_CellNoDis, MPI_Cell;

	// Create basic fields structure
	MPI_Type_create_struct(fields, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_CellNoDis);

	// Resize to cover alignment extent
	MPI_Aint lb, extent;
	MPI_Type_get_extent(MPI_CellNoDis, &lb, &extent);
	MPI_Type_create_resized(MPI_CellNoDis, lb, extent, &MPI_Cell);
	MPI_Type_commit(&MPI_Cell);

	/*Culture cells distribution
	if (columns > rows){
		int aux = rows;
		rows = columns;
		columns = aux;
	}*/

	if(nprocs > rows){
		nprocs = rows;
	}
	if(rank > nprocs){
		MPI_Finalize();
		return(0);
	}
	int num_culture_cells = rows * columns;

	/* Num of culture cells of every procesor */
	int self_rows = rows / nprocs;
	int offset_rows = rows - (self_rows * nprocs);

	/* Min and max culture cell for each procesor */
	int min_proc_rows = rank * self_rows + (rank < offset_rows ? rank : offset_rows);
	int max_proc_rows = min_proc_rows + self_rows + (rank < offset_rows ? 1 : 0);

	if (rank < offset_rows)
		self_rows++;

	/*Know neighbour processes */
	int upper_proc = (rank + 1) % nprocs;
	int lower_proc = rank - 1;
	if (lower_proc < 0)
		lower_proc = nprocs - 1;
	/*
	// Create the distributed surfaces for the process
	my_rows_number = rows / nprocs;
	rows_left = rows % nprocs;
	for (i = 0; i < rows_left; i++){
		if((i % nprocs) == rank)	rows_left;
	}

	my_begin;
	MPI_Status status;
	tag = 1000;
	if (nprocs > 1) {
		if (rank == 0){
			my_begin = 0;
			my_end = my_begin + my_rows_number - 1;
			MPI_Send( &my_end, 1, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
		}else{
			MPI_Recv( &my_begin, 1, MPI_INT, rank - 1, tag, MPI_COMM_WORLD, &status);
			my_end = my_begin + my_rows_number;
			if (rank != nprocs - 1){
				MPI_Send( &my_end, 1, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
			}
		}
	}
	*/
	printf("%d %d %d\n", rank, upper_proc, lower_proc);
	culture = (float *)malloc(sizeof(float) * (size_t)num_culture_cells);
	culture_cells = (short *)malloc(sizeof(short) * (size_t)num_culture_cells);
	short *culture_cells_aux = (short *)malloc(sizeof(short) * (size_t)num_culture_cells);
	float *culture_aux = (float *)malloc(sizeof(float) * (size_t)num_culture_cells);
	if (culture == NULL || culture_cells == NULL || culture_cells_aux == NULL || culture_aux == NULL)
	{
		fprintf(stderr, "-- Error allocating culture structures for size: %d x %d \n", rows, columns);
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}
#if !defined(CP_TABLON)
	timeInitCS = MPI_Wtime();
#endif
	for (i = 0; i < rows; i++)
		for (j = 0; j < columns; j++)
		{
			accessMat(culture, i, j) = 0.0;
			accessMat(culture_cells, i, j) = 0.0f;
		}

#if !defined(CP_TABLON)
	timeInitCS = MPI_Wtime() - timeInitCS;
#endif
#if !defined(CP_TABLON)
	timeInitCells = MPI_Wtime();
#endif
	/* Cell distribution */
	Cell *cells_subpart;
	if (rank == 0)
	{
		int num_init_cells_for_all[nprocs];
		for (i = 0; i < nprocs; i++)
		{
			num_init_cells_for_all[i] = 0;
		}
		int max_procs_rows[nprocs];
		MPI_Gather(&max_proc_rows, 1, MPI_INT, max_procs_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);

		for (i = 0; i < num_cells; i++)
		{
			cells[i].alive = true;
			// Initial age: Between 1 and 20
			cells[i].age = 1 + (int)(19 * erand48(cells[i].random_seq));
			// Initial storage: Between 10 and 20 units
			cells[i].storage = (float)(10 + 10 * erand48(cells[i].random_seq));
			// Initial position: Anywhere in the culture arena
			cells[i].pos_row = (float)(rows * erand48(cells[i].random_seq));
			cells[i].pos_col = (float)(columns * erand48(cells[i].random_seq));
			// Movement direction: Unity vector in a random direction
			cell_new_direction(&cells[i]);
			// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00
			cells[i].choose_mov[0] = 0.33f;
			cells[i].choose_mov[1] = 0.34f;
			cells[i].choose_mov[2] = 0.33f;
			j = 0;
			while (j < nprocs)
			{
				if ((int)cells[i].pos_row >= max_procs_rows[j])
				{
					j++;
				}
				else
				{

					num_init_cells_for_all[j]++;
					j = nprocs;
				}
			}
		}
		qsort(cells, num_cells, sizeof(Cell), &compareCells);
		MPI_Scatter(&num_init_cells_for_all, 1, MPI_INT, &num_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);
		cells_subpart = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
		if (cells_subpart == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}

		int displs[nprocs];
		int *sendcnts = num_init_cells_for_all;
		displs[0] = 0;
		for (i = 1; i < nprocs; i++)
		{
			displs[i] = displs[i - 1] + sendcnts[i - 1];
		}
		MPI_Scatterv(cells, sendcnts, displs, MPI_Cell, cells_subpart, num_cells, MPI_Cell, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Gather(&max_proc_rows, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(NULL, 1, MPI_INT, &num_cells, 1, MPI_INT, 0, MPI_COMM_WORLD);

		cells_subpart = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
		if (cells_subpart == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		MPI_Scatterv(0, 0, 0, MPI_Cell, cells_subpart, num_cells, MPI_Cell, 0, MPI_COMM_WORLD);
	}
	free(cells);
	cells = cells_subpart;
#if !defined(CP_TABLON)
	timeInitCells = MPI_Wtime() - timeInitCells;
#endif

	// Statistics: Initialize total number of cells, and max. alive
	sim_stat.history_total_cells = num_cells;
	sim_stat.history_max_alive_cells = num_cells;

#ifdef DEBUG
	/* Show initial cells data */
	printf("Initial cells data: %d\n", num_cells);
	for (i = 0; i < num_cells; i++)
	{
		printf("\tCell %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
			   i,
			   cells[i].pos_row,
			   cells[i].pos_col,
			   cells[i].mov_row,
			   cells[i].mov_col,
			   cells[i].choose_mov[0],
			   cells[i].choose_mov[1],
			   cells[i].choose_mov[2],
			   cells[i].storage,
			   cells[i].age);
	}
#endif // DEBUG

/* 4. Simulation */
#if !defined(CP_TABLON)
	timeML = MPI_Wtime();
#endif
	float current_max_food = 0.0f;
	int num_cells_alive = num_cells;
	int iter;
	for (iter = 0; iter < max_iter && current_max_food <= max_food && any_cell_alive; iter++)
	{
		if (rank == 0)
			printf("%d\n", iter);

		int step_new_cells = 0;
		int step_dead_cells = 0;

		/* 4.1. Spreading new food */
		// Across the whole culture
		int num_new_sources = (int)(rows * columns * food_density);
#if !defined(CP_TABLON)
		timeNormalSpreadingL = MPI_Wtime();
#endif
		for (i = 0; i < num_new_sources; i++)
		{
			int row = (int)(rows * erand48(food_random_seq));
			int col = (int)(columns * erand48(food_random_seq));
			float food = (float)(food_level * erand48(food_random_seq));
			accessMat(culture, row, col) = accessMat(culture, row, col) + food;
		}
#if !defined(CP_TABLON)
		timeNormalSpreadingL = MPI_Wtime() - timeNormalSpreadingL;
		timeNormalSpreadingT += timeNormalSpreadingL;
#endif
#if !defined(CP_TABLON)
		timeSpecialSpreadingL = MPI_Wtime();
#endif
		// In the special food spot
		if (food_spot_active)
		{
			num_new_sources = (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density);
			for (i = 0; i < num_new_sources; i++)
			{
				int row = food_spot_row + (int)(food_spot_size_rows * erand48(food_spot_random_seq));
				int col = food_spot_col + (int)(food_spot_size_cols * erand48(food_spot_random_seq));
				float food = (float)(food_spot_level * erand48(food_spot_random_seq));
				accessMat(culture, row, col) = accessMat(culture, row, col) + food;
			}
		}
#if !defined(CP_TABLON)
		timeSpecialSpreadingL = MPI_Wtime() - timeSpecialSpreadingL;
		timeSpecialSpreadingT += timeSpecialSpreadingL;
#endif

		/* 4.2. Prepare ancillary data structures */
		/* 4.2.1. Clear ancillary structure of the culture to account alive cells in a position after movement */

		/* 4.2.2. Allocate ancillary structure to store the food level to be shared by cells in the same culture place */
		float *food_to_share = (float *)malloc(sizeof(float) * num_cells);
		if (culture == NULL || culture_cells == NULL)
		{
			fprintf(stderr, "-- Error allocating culture structures for size: %d x %d \n", rows, columns);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}

/* 4.3. Cell movements */
#if !defined(CP_TABLON)
		timeCellMovementL = MPI_Wtime();
#endif
		int total_to_upper = 0;
		int total_to_lower = 0;
		Cell *cells_to_upper = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
		Cell *cells_to_lower = (Cell *)malloc(sizeof(Cell) * (size_t)num_cells);
		if (cells_to_upper == NULL || cells_to_lower == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		for (i = 0; i < num_cells; i++)
		{
			bool to_upper = false, to_lower = false;
			cells[i].age++;
			// Statistics: Max age of a cell in the simulation history
			if (cells[i].age > sim_stat.history_max_age)
				sim_stat.history_max_age = cells[i].age;

			/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
			if (cells[i].storage < 0.1f)
			{
				// Cell has died
				cells[i].alive = false;
				num_cells_alive--;
				step_dead_cells++;
				continue;
			}
			if (cells[i].storage < 1.0f)
			{
				// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
				cells[i].storage -= 0.2f;
			}
			else
			{
				// Consume energy to move
				cells[i].storage -= 1.0f;

				/* 4.3.2. Choose movement direction */
				float prob = (float)erand48(cells[i].random_seq);
				if (prob < cells[i].choose_mov[0])
				{
					// Turn left (90 degrees)
					float tmp = cells[i].mov_col;
					cells[i].mov_col = cells[i].mov_row;
					cells[i].mov_row = -tmp;
				}
				else if (prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1])
				{
					// Turn right (90 degrees)
					float tmp = cells[i].mov_row;
					cells[i].mov_row = cells[i].mov_col;
					cells[i].mov_col = -tmp;
				}
				// else do not change the direction

				/* 4.3.3. Update position moving in the choosen direction*/
				cells[i].pos_row += cells[i].mov_row;
				cells[i].pos_col += cells[i].mov_col;

				if (cells[i].pos_row < min_proc_rows)
					to_lower = true;
				else if (cells[i].pos_row >= max_proc_rows){
					to_upper = true;
				}
				//printf("rank %d %lf %lf %d %d \n",rank, cells[i].pos_row,cells[i].pos_col,(int)cells[i].pos_row,max_proc_rows);

				// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
				if (cells[i].pos_row < 0)
					cells[i].pos_row += rows;
				if (cells[i].pos_row >= rows)
					cells[i].pos_row -= rows;
				if (cells[i].pos_col < 0)
					cells[i].pos_col += columns;
				if (cells[i].pos_col >= columns)
					cells[i].pos_col -= columns;


			}

			if (to_upper)
			{
				cells_to_upper[total_to_upper] = cells[i];
				cells[i].alive = false;

				total_to_upper++;
				num_cells_alive--;
			}
			else if (to_lower)
			{
				cells_to_lower[total_to_lower] = cells[i];
				cells[i].alive = false;

				total_to_lower++;
				num_cells_alive--;
			}
			else
			{
				/* 4.3.4. Annotate that there is one more cell in this culture position */
				accessMat(culture_cells, cells[i].pos_row, cells[i].pos_col) += 1;
				/* 4.3.5. Annotate the amount of food to be shared in this culture position */
				food_to_share[i] = accessMat(culture, cells[i].pos_row, cells[i].pos_col);
			}
		} // End cell movements
		/* Sync the culture_cells in all processes */
		MPI_Allreduce(culture_cells, culture_cells_aux, num_culture_cells, MPI_SHORT, MPI_SUM, MPI_COMM_WORLD);
		short *aux_for_interchange_culture_cells = culture_cells;
		culture_cells = culture_cells_aux;
		culture_cells_aux = aux_for_interchange_culture_cells;

		printf("ini %d iter %d\n", rank, iter);
		MPI_Request req;
		MPI_Status stat;
		int size_from_upper, size_from_lower;
		MPI_Isend(cells_to_upper, total_to_upper, MPI_Cell, upper_proc, 0, MPI_COMM_WORLD, &req);
		for (i = 0; i < total_to_upper; i++)
		{
			printf("\tsend iter %d to upper Rank %d Cell %d of %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				   iter, rank, i, total_to_upper,
				   cells_to_upper[i].pos_row,
				   cells_to_upper[i].pos_col,
				   cells_to_upper[i].mov_row,
				   cells_to_upper[i].mov_col,
				   cells_to_upper[i].choose_mov[0],
				   cells_to_upper[i].choose_mov[1],
				   cells_to_upper[i].choose_mov[2],
				   cells_to_upper[i].storage,
				   cells_to_upper[i].age);
		}
		MPI_Probe(lower_proc, 0, MPI_COMM_WORLD, &stat);
		MPI_Get_count(&stat, MPI_Cell, &size_from_lower);
		free(cells_to_upper);
		Cell *cells_from_lower = (Cell *)malloc(sizeof(Cell) * (size_t)size_from_lower);
		if (cells_from_lower == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		printf("pre recv from lower %d iter %d\n", rank, iter);

		MPI_Recv(cells_from_lower, size_from_lower, MPI_Cell, lower_proc, 0, MPI_COMM_WORLD, &stat);
		printf("post recv from lower %d iter %d\n", rank, iter);

		for (i = 0; i < size_from_lower; i++)
		{
			printf("\titer %d From lower Rank %d Cell %d of %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				   iter, rank, i, size_from_lower,
				   cells_from_lower[i].pos_row,
				   cells_from_lower[i].pos_col,
				   cells_from_lower[i].mov_row,
				   cells_from_lower[i].mov_col,
				   cells_from_lower[i].choose_mov[0],
				   cells_from_lower[i].choose_mov[1],
				   cells_from_lower[i].choose_mov[2],
				   cells_from_lower[i].storage,
				   cells_from_lower[i].age);
		}

		//MPI_Wait(&req, &stat);
				printf( "pre send to lower %d iter %d\n", rank, iter);

		MPI_Isend(cells_to_lower, total_to_lower, MPI_Cell, lower_proc, 0, MPI_COMM_WORLD, &req);

		for (i = 0; i < total_to_lower; i++)
		{
			printf("\tsend iter %d to lower Rank %d Cell %d of %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				   iter, rank, i, total_to_lower,
				   cells_to_lower[i].pos_row,
				   cells_to_lower[i].pos_col,
				   cells_to_lower[i].mov_row,
				   cells_to_lower[i].mov_col,
				   cells_to_lower[i].choose_mov[0],
				   cells_to_lower[i].choose_mov[1],
				   cells_to_lower[i].choose_mov[2],
				   cells_to_lower[i].storage,
				   cells_to_lower[i].age);
		}

		MPI_Probe(upper_proc, 0, MPI_COMM_WORLD, &stat);
		MPI_Get_count(&stat, MPI_Cell, &size_from_upper);
		for (i = 0; i < total_to_lower; i++)
			cells_to_lower[i].alive = false;
		free(cells_to_lower);
		Cell *cells_from_upper = (Cell *)malloc(sizeof(Cell) * (size_t)size_from_upper);
		if (cells_from_upper == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}

				printf("pre recv from upper %d iter %d and size_from_upper %d num_cells %d\n", rank, iter, size_from_upper, num_cells);

		MPI_Recv(cells_from_upper, size_from_upper, MPI_Cell, upper_proc, 0, MPI_COMM_WORLD, &stat);
						printf("post recv from upper %d iter %d\n", rank, iter);
	//	MPI_Wait(&req,  &stat);

		for (i = 0; i < size_from_upper; i++)
		{
			printf("\t iter %d From upper Rank %d Cell %d of %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				   iter, rank, i, size_from_upper,
				   cells_from_upper[i].pos_row,
				   cells_from_upper[i].pos_col,
				   cells_from_upper[i].mov_row,
				   cells_from_upper[i].mov_col,
				   cells_from_upper[i].choose_mov[0],
				   cells_from_upper[i].choose_mov[1],
				   cells_from_upper[i].choose_mov[2],
				   cells_from_upper[i].storage,
				   cells_from_upper[i].age);
		}
			//	printf("rank %d antiguo num cells %d\n", rank, num_cells);
		cells = (Cell *)realloc(cells, sizeof(Cell) * (num_cells + size_from_lower + size_from_upper));
		food_to_share = (float *)realloc(food_to_share, sizeof(float) * (num_cells + size_from_lower + size_from_upper));

		if (cells == NULL || food_to_share == NULL)
		{
			fprintf(stderr, "-- Error allocating: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
				printf("pre order from upper %d iter %d\n", rank, iter);

		for (i = 0; i < size_from_upper; i++)
		{
			if((int)cells_from_upper[i].pos_row < min_proc_rows || (int)cells_from_upper[i].pos_row >= max_proc_rows)
				printf("CELULA MAL MANDADA FROM UPPER iter %d rank %d cell %d %d < %d < %d\n", iter, rank, i, min_proc_rows,(int)cells_from_upper[i].pos_row,max_proc_rows);
			cells[num_cells + i] = cells_from_upper[i];
			cells[num_cells + i].alive = true;

			num_cells_alive++;
			/* 4.3.4. Annotate that there is one more cell in this culture position */
			accessMat(culture_cells, cells_from_upper[i].pos_row, cells_from_upper[i].pos_col) += 1;
			/* 4.3.5. Annotate the amount of food to be shared in this culture position */
			food_to_share[num_cells + i] = accessMat(culture, cells_from_upper[i].pos_row, cells_from_upper[i].pos_col);

		}
		free(cells_from_upper);
						printf("pre order from lower %d iter %d\n", rank, iter);

		for (i = 0; i < size_from_lower; i++)
		{
			if((int)cells_from_lower[i].pos_row < min_proc_rows || (int)cells_from_lower[i].pos_row >= max_proc_rows)
				printf("CELULA MAL MANDADA FROM LOWER iter %d rank %d cell %d %d < %d < %d\n", iter, rank, i, min_proc_rows,(int)cells_from_lower[i].pos_row,max_proc_rows);
			cells[num_cells + size_from_upper + i] = cells_from_lower[i];
			cells[num_cells + size_from_upper + i].alive = true;

			num_cells_alive++;
			/* 4.3.4. Annotate that there is one more cell in this culture position */
			accessMat(culture_cells, cells_from_lower[i].pos_row, cells_from_lower[i].pos_col) += 1;
			/* 4.3.5. Annotate the amount of food to be shared in this culture position */
			food_to_share[num_cells + size_from_upper + i] = accessMat(culture, cells_from_lower[i].pos_row, cells_from_lower[i].pos_col);

		}
		free(cells_from_lower);
		printf("fin %d iter %d\n", rank, iter);

		num_cells = num_cells + size_from_lower + size_from_upper;

#if !defined(CP_TABLON)
		timeCellMovementL = MPI_Wtime() - timeCellMovementL;
		timeCellMovementT += timeCellMovementL;
#endif
/* 4.4. Cell actions */
#if !defined(CP_TABLON)
		timeCellActionsL = MPI_Wtime();
#endif
		// Space for the list of new cells (maximum number of new cells is num_cells)
		Cell *new_cells = (Cell *)malloc(sizeof(Cell) * num_cells);
		if (new_cells == NULL)
		{
			fprintf(stderr, "-- Error allocating new cells structures for: %d cells\n", num_cells);
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}

		for (i = 0; i < num_cells; i++)
		{
			if (cells[i].alive)
			{
				/* 4.4.1. Food harvesting */
				float food = food_to_share[i];
				short count = accessMat(culture_cells, cells[i].pos_row, cells[i].pos_col);
				float my_food = food / count;
				cells[i].storage += my_food;

				/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
				if (cells[i].age > 30 && cells[i].storage > 20)
				{
					// Split: Create new cell
					num_cells_alive++;
					sim_stat.history_total_cells++;
					step_new_cells++;

					// New cell is a copy of parent cell
					new_cells[step_new_cells - 1] = cells[i];

					// Split energy stored and update age in both cells
					cells[i].storage /= 2.0f;
					new_cells[step_new_cells - 1].storage /= 2.0f;
					cells[i].age = 1;
					new_cells[step_new_cells - 1].age = 1;

					// Random seed for the new cell, obtained using the parent random sequence
					new_cells[step_new_cells - 1].random_seq[0] = (unsigned short)nrand48(cells[i].random_seq);
					new_cells[step_new_cells - 1].random_seq[1] = (unsigned short)nrand48(cells[i].random_seq);
					new_cells[step_new_cells - 1].random_seq[2] = (unsigned short)nrand48(cells[i].random_seq);

					// Both cells start in random directions
					cell_new_direction(&cells[i]);
					cell_new_direction(&new_cells[step_new_cells - 1]);

					// Mutations of the movement genes in both cells
					cell_mutation(&cells[i]);
					cell_mutation(&new_cells[step_new_cells - 1]);
				}
				accessMat(culture, cells[i].pos_row, cells[i].pos_col) = 0.0f;
			}
		} // End cell actions
#if !defined(CP_TABLON)
		timeCellActionsL = MPI_Wtime() - timeCellActionsL;
		timeCellActionsT += timeCellActionsL;
#endif

/* 4.6. Clean dead cells from the original list */
// 4.6.1. Move alive cells to the left to substitute dead cells
#if !defined(CP_TABLON)
		timeMovingAliveCellsL = MPI_Wtime();
#endif
		int free_position = 0;
		int alive_in_main_list = 0;
		for (i = 0; i < num_cells; i++)
		{
			if (cells[i].alive)
			{
				alive_in_main_list++;
				if (free_position != i)
				{
					cells[free_position] = cells[i];
				}
				free_position++;
			}
		}

#if !defined(CP_TABLON)
		timeMovingAliveCellsL = MPI_Wtime() - timeMovingAliveCellsL;
		timeMovingAliveCellsT += timeMovingAliveCellsL;
#endif
		// 4.6.2. Reduce the storage space of the list to the current number of cells
		num_cells = alive_in_main_list;
		cells = (Cell *)realloc(cells, sizeof(Cell) * num_cells);

/* 4.7. Join cell lists: Old and new cells list */
#if !defined(CP_TABLON)
		timeJoinCellsListL = MPI_Wtime();
#endif
		if (step_new_cells > 0)
		{
			cells = (Cell *)realloc(cells, sizeof(Cell) * (num_cells + step_new_cells));
			for (j = 0; j < step_new_cells; j++)
				cells[num_cells + j] = new_cells[j];
			num_cells += step_new_cells;
		}
#if !defined(CP_TABLON)
		timeJoinCellsListL = MPI_Wtime() - timeJoinCellsListL;
		timeJoinCellsListT += timeJoinCellsListL;
#endif
		free(new_cells);

/* 4.8. Decrease non-harvested food */
#if !defined(CP_TABLON)
		timeDecreaseFoodL = MPI_Wtime();
#endif
		/* Sync the culture in all processes */
		MPI_Allreduce(culture, culture_aux, num_culture_cells, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		float *aux_for_interchange_culture = culture;
		culture = culture_aux;
		culture_aux = aux_for_interchange_culture;

		current_max_food = 0.0f;
		for (i = 0; i < rows; i++)
			for (j = 0; j < columns; j++)
			{
				accessMat(culture_cells, i, j) = 0.0f;

				accessMat(culture, i, j) *= 0.95f; // Reduce 5%
				if (accessMat(culture, i, j) > current_max_food)
					current_max_food = accessMat(culture, i, j);
			}

#if !defined(CP_TABLON)
		timeDecreaseFoodL = MPI_Wtime() - timeDecreaseFoodL;
		timeDecreaseFoodT += timeDecreaseFoodL;
#endif

#if !defined(CP_TABLON)
		timeDataCollectionL = MPI_Wtime();
#endif

		int aux = 0;
		/*Check if any cell is still alive */
		if (rank == 0)
		{
			MPI_Reduce(&num_cells_alive, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			if (aux == 0)
			{
				any_cell_alive = false;
			}
			MPI_Bcast(&any_cell_alive, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Reduce(&num_cells_alive, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&any_cell_alive, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
		}
		/* 4.9. Statistics */
		// Statistics: Max food
		float aux_float;
		MPI_Reduce(&current_max_food, &aux_float, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		if (aux_float > sim_stat.history_max_food)
			sim_stat.history_max_food = aux_float;
		// Statistics: Max new cells per step
		MPI_Reduce(&step_new_cells, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (aux > sim_stat.history_max_new_cells)
			sim_stat.history_max_new_cells = aux;
		// Statistics: Accumulated dead and Max dead cells per step
		sim_stat.history_dead_cells += step_dead_cells;
		MPI_Reduce(&step_dead_cells, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (aux > sim_stat.history_max_dead_cells)
			sim_stat.history_max_dead_cells = aux;
		// Statistics: Max alive cells per step
		MPI_Reduce(&num_cells_alive, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (aux > sim_stat.history_max_alive_cells)
			sim_stat.history_max_alive_cells = aux;

#if !defined(CP_TABLON)
		timeDataCollectionL = MPI_Wtime() - timeDataCollectionL;
		timeDataCollectionT += timeDataCollectionL;
#endif
#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		//if (rank == 0)

			//print_status(iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat);
#endif // DEBUG
		}
#if !defined(CP_TABLON)
		timeML = MPI_Wtime() - timeML;
#endif

		int aux = 0;
		MPI_Reduce(&sim_stat.history_dead_cells, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		sim_stat.history_dead_cells = aux;
		MPI_Reduce(&sim_stat.history_max_age, &aux, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
		sim_stat.history_max_age = aux;
		MPI_Reduce(&sim_stat.history_total_cells, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		sim_stat.history_total_cells = aux;

		MPI_Reduce(&num_cells_alive, &aux, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		num_cells_alive = aux;
		/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

		/* 5. Stop global time */
		MPI_Barrier(MPI_COMM_WORLD);
		ttotal = cp_Wtime() - ttotal;

#ifdef DEBUG
		printf("List of cells at the end of the simulation: %d\n\n", num_cells);
		for (i = 0; i < num_cells; i++)
		{
			printf("Cell %d, Alive: %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				   i,
				   cells[i].alive,
				   cells[i].pos_row,
				   cells[i].pos_col,
				   cells[i].mov_row,
				   cells[i].mov_col,
				   cells[i].choose_mov[0],
				   cells[i].choose_mov[1],
				   cells[i].choose_mov[2],
				   cells[i].storage,
				   cells[i].age);
		}
#endif // DEBUG
#if !defined(CP_TABLON)
		if (rank == 0)
		{
			// 6.1.1 Disgragated time used on each loop
			double timeGatherer[nprocs];
			double totalGatheredTime;
			int i;
			printf("---------TIMES---------\n");

			MPI_Gather(&timeInitCS, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Init culture:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeInitCells, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Init cells:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeML, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Main Loop:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeNormalSpreadingT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Normal food spreading:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeSpecialSpreadingT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Special food spreading:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeCellMovementT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Cell movement:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeCellActionsT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Cell action:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeMovingAliveCellsT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Reposition alive cells:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeJoinCellsListT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Join cells:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeDecreaseFoodT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Decrease food:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);

			MPI_Gather(&timeDataCollectionT, 1, MPI_DOUBLE, timeGatherer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			printf("\n" ANSI_COLOR_RED "Statistics Sync:\n\t" ANSI_COLOR_RESET);
			printTimes(timeGatherer, nprocs, ttotal);
		}
		else
		{
			MPI_Gather(&timeInitCS, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeInitCells, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeML, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeNormalSpreadingT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeSpecialSpreadingT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeCellMovementT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeCellActionsT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeMovingAliveCellsT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeJoinCellsListT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeDecreaseFoodT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gather(&timeDataCollectionT, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
#endif

		/* 6. Output for leaderboard */
		if (rank == 0)
		{
			printf("\n");
			/* 6.1. Total computation time */
			printf("Time: %lf\n", ttotal);

			/* 6.2. Results: Number of iterations and other statistics */
			printf("Result: %d, ", iter);
			printf("%d, %d, %d, %d, %d, %d, %d, %f\n",
				   num_cells_alive,
				   sim_stat.history_total_cells,
				   sim_stat.history_dead_cells,
				   sim_stat.history_max_alive_cells,
				   sim_stat.history_max_new_cells,
				   sim_stat.history_max_dead_cells,
				   sim_stat.history_max_age,
				   sim_stat.history_max_food);
		}
		/* 7. Free resources */
		free(culture);
		free(culture_cells);
		free(cells);

		/* 8. End */
		MPI_Finalize();
		return 0;
	}
	/* function to compare cells by position*/
	int compareCells(const void *_a, const void *_b)
	{

		Cell *a, *b;

		a = (Cell *)_a;
		b = (Cell *)_b;
		if (a->pos_row < b->pos_row)
		{
			return -1;
		}
		if (a->pos_row > b->pos_row)
		{
			return 1;
		}
		if (a->pos_col < b->pos_col)
		{
			return -1;
		}
		return 1;
	}
