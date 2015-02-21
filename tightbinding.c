#include <tb2.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "doublecomplex.h"
#include "hamiltonian.h"
#include "utils.h"
#include "constants.h"
#include "tightbinding.h"
#include "plane.h"
#include "lattice.h"

const double TB_MN_X = 0.042;

static tb_geometry current_tb_geometry = -1;
static tb_parameterization current_tb_parameterization = -1;
static tb_material current_tb_material = -1;

const char* tb2_input_filename_prefix = "werdi_input_tb_";

static int get_bulk_dimension(tb_geometry g)
{
	switch(g)
	{
		case GEOM_XY:
			return 40;
			break;
		case GEOM_Z:
			return 80;
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", g);
			abort();
	}
}

static int get_bulk_p_bands_cnt(tb_geometry g)
{
	switch(g)
	{
		case GEOM_XY:
			return 6;
			break;
		case GEOM_Z:
			return 12;
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", g);
			abort();
	}
}

static int get_bulk_p_bands_start(tb_parameterization p, tb_geometry g)
{
	switch(g)
	{
		case GEOM_XY:
			return 2;
			break;
		case GEOM_Z:
			return 4;
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", g);
			abort();
	}
}

static int get_periodic_bc_dim(tb_material m)
{
	switch(m)
	{
		case BULK:
			return 3;
			break;
		case PLANE:
			return 2;
			break;
		default:
			fprintf(stderr, "Unsupported material: %d\n", m);
			abort();
	}
}

static const char* generate_tb2_input_filename(tb_parameterization p, tb_geometry g, tb_material m)
{
	char* filename = alloc_memory(sizeof(char) * 256);
	strcpy(filename, tb2_input_filename_prefix);
	switch(p)
	{
		case JANCU:
			strcat(filename, "jancu");
			break;
		case DICARLO:
			strcat(filename, "dicarlo");
			break;
		default:
			fprintf(stderr, "Unsuported parameterization: %d\n", p);
			abort();
	}
	strcat(filename, "_");
	switch(g)
	{
		case GEOM_XY:
			strcat(filename, "xy");
			break;
		case GEOM_Z:
			strcat(filename, "z");
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", g);
			abort();
	}
	strcat(filename, "_");
	switch(m)
	{
		case BULK:
			strcat(filename, "bulk");
			break;
		case PLANE:
			strcat(filename, "plane");
			break;
		default:
			fprintf(stderr, "Unsupported material: %d\n", m);
			break;
	}
	strcat(filename, ".dat");
	printf("%s\n", filename);
	fflush(stdout);
	return filename;
}

const int STARTING_PLANE = 3; //3. kolumna inputu, program robi rozne rzeczy z 1. warstwa GaMnAs, wiec zaczynamy od kolejnej (2. w GaMnAs, 3. w ogole)
const double BASE_POTENTIAL = 0; // electric potential in x = (0,0,0)

/* K-step for tight binding, in TB units */
double find_cstep_tb(ham_struct* this, double pmax, int npts)
{
        return 1.0/(double)npts;
}
/* Integrands for tight binding, in 3D and 2D */
double bril_zone_integrand_tb_3d(ham_struct* this, double (*int_fun) (double, double, double), const double* x)
{
        return int_fun(x[0] - 0.5, x[1] - 0.5, x[2] - 0.5);
}
double bril_zone_integrand_tb_2d(ham_struct* this, double (*int_fun) (double, double, double), const double* x)
{
        return int_fun(x[0] - 0.5, x[1] - 0.5, 0);
}

static double electric_potential(double V0, double pos[3], double electric_field[3])
{
    return V0 + electric_field[0]*pos[0] + electric_field[1]*pos[1] + electric_field[2]*pos[2];
}

void calculate_tb2_bulk_hamiltonian(ham_struct * this, double c1, double c2, double c3, doublecomplex matrix[][])
{
	double c[3];
	c[0] = c1;
	c[1] = c2;
	c[2] = c3;
	tb2setmagdir(this->magdir);
	tb2hamilton(this->D / Delta0(this), c, matrix); // D/Delta0 = magfac
}

void calculate_tb2_plane_hamiltonian(ham_struct* this, double c1, double c2, double c3, doublecomplex *matrix)
{
    double c[3];
    c[0] = c1;
    c[1] = c2;
    c[2] = c3;
    
    double delta_factor = this->D / Delta0(this);
    plane_ham_struct* phs = (ham_struct*) this;
    int nlayers = phs->nlayers; // Number of planes (layers)
    doublecomplex* work = phs->single_plane_matrix;
    int dim = this->dim; // Dimension of multi-plane Hamiltonian
    int planehamdim = dim / nlayers; // Dimension of the single-plane Hamiltonian
    tb2setmagdir(this->magdir);
    
    for (int p = 0; p < nlayers; p++)
    {
        tb2hamiltonplane(delta_factor, c, p + STARTING_PLANE, 0, 0, work);
        tb2atompositions(p + STARTING_PLANE, phs->atom_positions);
        for (int k = 0; k < planehamdim; k++)
        {
            for (int l = 0; l < planehamdim; l++)
            {
                matrix[(p * planehamdim + k) * dim + p * planehamdim + l] = work[k * planehamdim + l];
            }
            matrix[(p * planehamdim + k) * dim + p * planehamdim + k].real += electric_potential(BASE_POTENTIAL, phs->atom_positions + 3*k, phs->electric_field) //
                + (phs->concentration_gradient == NULL ? 0 : phs->concentration_gradient(phs->atom_positions + 3 * k));
        }
        if (p < nlayers - 1) 
        {
            tb2hamiltonplane(delta_factor, c, p + STARTING_PLANE, 1, 0, work);
            for (int k = 0; k < planehamdim; k++)
            {
                for (int l = 0; l < planehamdim; l++)
                {
                    int source_idx = k * planehamdim + l;
                    matrix[((p + 1) * planehamdim + k) * dim + p * planehamdim + l] = work[source_idx];
                    int transv_dest_idx = (p * planehamdim + l) * dim + (p + 1) * planehamdim + k;
                    matrix[transv_dest_idx].real = work[source_idx].real;
                    matrix[transv_dest_idx].imag = - work[source_idx].imag;
                }
            }
        }        
        for (int q = p + 2; q < nlayers; q++)
        {
            for (int k = 0; k < planehamdim; k++) 
            {
                for (int l = 0; l < planehamdim; l++) 
                {
                    matrix[(q * planehamdim + k) * dim + p * planehamdim + l] = CMPLX_ZERO;
                    matrix[(p * planehamdim + k) * dim + q * planehamdim + l] = CMPLX_ZERO;
                }
            }
        }
    }
}

void exact_tb_plane_derivative(ham_struct* this, double c1, double c2, double c3, double dir1, double dir2, double dir3, doublecomplex *dHdc)
{
    double c[3];
    c[0] = c1;
    c[1] = c2;
    c[2] = c3;
    double dir[3];
    dir[0] = dir1;
    dir[1] = dir2;
    dir[2] = dir3;
    
    double delta_factor = this->D / Delta0(this);
    plane_ham_struct* phs = (ham_struct*) this;
    int nlayers = phs->nlayers; // Number of planes (layers)
    doublecomplex* work = phs->single_plane_matrix;
    int dim = this->dim; // Dimension of multi-plane Hamiltonian
    int planehamdim = dim / nlayers; // Dimension of the single-plane Hamiltonian
    
    for (int p = 0; p < nlayers; p++)
    {
        tb2hamiltonplane_derivative(delta_factor, c, dir, p + STARTING_PLANE, 0, work);
        for (int k = 0; k < planehamdim; k++)
        {
            for (int l = 0; l < planehamdim; l++)
            {
                dHdc[(p * planehamdim + k) * dim + p * planehamdim + l] = work[k * planehamdim + l];
            }
        }
        if (p < nlayers - 1) 
        {
            tb2hamiltonplane_derivative(delta_factor, c, dir, p + STARTING_PLANE, 1, work);
            for (int k = 0; k < planehamdim; k++)
            {
                for (int l = 0; l < planehamdim; l++)
                {
                    int source_idx = k * planehamdim + l;
                    dHdc[((p + 1) * planehamdim + k) * dim + p * planehamdim + l] = work[source_idx];
                    int transv_dest_idx = (p * planehamdim + l) * dim + (p + 1) * planehamdim + k;
                    dHdc[transv_dest_idx].real = work[source_idx].real;
                    dHdc[transv_dest_idx].imag = - work[source_idx].imag;
                }
            }
        }        
        for (int q = p + 2; q < nlayers; q++)
        {
            for (int k = 0; k < planehamdim; k++) 
            {
                for (int l = 0; l < planehamdim; l++) 
                {
                    dHdc[(q * planehamdim + k) * dim + p * planehamdim + l] = CMPLX_ZERO;
                    dHdc[(p * planehamdim + k) * dim + q * planehamdim + l] = CMPLX_ZERO;
                }
            }
        }
    }
}


static doublecomplex* create_S_z_bulk(int dim)
{
	doublecomplex* sz = initialize_cmplx_matrix(dim);	
	int half_dim = dim / 2;
	for (int i = 0; i < half_dim; i++)
	{
		sz[i * dim + i].real = 0.5;
		sz[(i + half_dim) * dim + i + half_dim].real = - 0.5;
	}

	return sz;
}

static doublecomplex* create_S_x_bulk(int dim)
{
	doublecomplex* sx = initialize_cmplx_matrix(dim);	
	int half_dim = dim / 2;
	for (int i = 0; i < half_dim; i++)
	{
		sx[(i + half_dim) * dim + i].real = 0.5;
		sx[i * dim + i + half_dim].real = 0.5;
	}
	return sx;
}

static doublecomplex* create_S_y_bulk(int dim)
{
	doublecomplex* sy = initialize_cmplx_matrix(dim);	
	int half_dim = dim / 2;
	for (int i = 0; i < half_dim; i++)
	{
		sy[i * dim + i + half_dim].imag = 0.5;
		sy[(i + half_dim) * dim + i].imag = - 0.5;
	}
	return sy;
}

static doublecomplex* create_S_minus(const doublecomplex* sx, const doublecomplex* sy, int dim)
{
	doublecomplex* sminus = initialize_cmplx_matrix(dim);
	double sqrt2 = sqrt(2.0);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			int idx = i * dim + j;
			sminus[idx].real = (sx[idx].real + sy[idx].imag) / sqrt2;
			sminus[idx].imag = (sx[idx].imag - sy[idx].real) / sqrt2;
		}
	}
	return sminus;
}

static doublecomplex* create_S_plus(const doublecomplex* sx, const doublecomplex* sy, int dim)
{
	doublecomplex* splus = initialize_cmplx_matrix(dim);
	double sqrt2 = sqrt(2.0);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			int idx = i * dim + j;
			splus[idx].real = (sx[idx].real - sy[idx].imag) / sqrt2;
			splus[idx].imag = (sx[idx].imag + sy[idx].real) / sqrt2;
		}
	}
	return splus;
}


static doublecomplex* create_plane_spin_matrix_from_bulk(doublecomplex* bulk_sm, int bulk_dim, int nlayers)
{
	int dim = bulk_dim * nlayers;
	doublecomplex* sm = initialize_cmplx_matrix(dim);

	for (int l = 0; l < nlayers; l++)
	{
		for (int i = 0; i < bulk_dim; i++)
		{
			int i2 = l * bulk_dim + i;
			for (int j = 0; j < bulk_dim; j++)
			{
				int j2 = l * bulk_dim + j;
				sm[i2 * dim + j2] = bulk_sm[i * bulk_dim + j];
			}
		}
	}
	
	return sm;
}

/*
   Create a S_z matrix for TB plane hamiltonian.
   bulk_dim - dimension of corresponding bulk hamiltonian
   nlayers - number of layers
 */
static doublecomplex* create_S_z_plane(int bulk_dim, int nlayers)
{
	if (bulk_dim % 2 != 0)
	{
		fprintf(stderr, "Odd dimension for TB hamiltonian");
		abort();
	}
	doublecomplex* bulk_sz = create_S_z_bulk(bulk_dim);
	doublecomplex* sz = create_plane_spin_matrix_from_bulk(bulk_sz, bulk_dim, nlayers);
	free(bulk_sz);
	return sz;
}


/*
   Create a S_x matrix for TB plane hamiltonian.
   bulk_dim - dimension of corresponding bulk hamiltonian
   nlayers - number of layers
 */
static doublecomplex* create_S_x_plane(int bulk_dim, int nlayers)
{
	if (bulk_dim % 2 != 0)
	{
		fprintf(stderr, "Odd dimension for TB hamiltonian");
		abort();
	}
	doublecomplex* bulk_sx = create_S_x_bulk(bulk_dim);
	doublecomplex* sx = create_plane_spin_matrix_from_bulk(bulk_sx, bulk_dim, nlayers);
	free(bulk_sx);
	return sx;
}

/*
   Create a S_y matrix for TB plane hamiltonian.
   bulk_dim - dimension of corresponding bulk hamiltonian
   nlayers - number of layers
 */
static doublecomplex* create_S_y_plane(int bulk_dim, int nlayers)
{
	if (bulk_dim % 2 != 0)
	{
		fprintf(stderr, "Odd dimension for TB hamiltonian");
		abort();
	}
	doublecomplex* bulk_sy = create_S_y_bulk(bulk_dim);
	doublecomplex* sy = create_plane_spin_matrix_from_bulk(bulk_sy, bulk_dim, nlayers);
	free(bulk_sy);
	return sy;
}


void destroy_tb_plane_hamiltonian_structure(ham_struct *this)
{
	plane_ham_struct* phs = (plane_ham_struct*) this;
	destroy_plane_ham_struct(phs);
	destroy_ham_struct(&phs->super);
}

static void setup_xy_geometry(ham_struct* hs)
{
	const double aSI[3][3] = { { 0.28165, 0.28165, 0 }, { 0.28165, 0, 0.28165 }, { 0, 0.28165, 0.28165 } };

	hs->L = 0.5633 / 2 / M_PI / sqrt(3) ;  // Length of b1 for zero strain.

	setup_geometry(hs, aSI);

	/*
	// Values for zero strain (exx == 0)
	//
	
       double aa = sqrt(3) / 2;
       hs->kx_vect[0] = aa;
       hs->kx_vect[1] = aa;
       hs->kx_vect[2] = 0;
       hs->ky_vect[0] = aa;
       hs->ky_vect[1] = 0;
       hs->ky_vect[2] = aa;
       hs->kz_vect[0] = 0;
       hs->kz_vect[1] = aa;
       hs->kz_vect[2] = aa;

       hs->bvect_lengths[0] = 1;
       hs->bvect_lengths[1] = 1;
       hs->bvect_lengths[2] = 1;
	
       hs->volume_factor = 1 / sqrt(2);	
	hs->volume_factor_12 = sqrt(3) / 2;
	hs->cosine_product = 2 * sqrt(2) / 3 / sqrt(3);
	hs->cosine_product_12 = 2.0 / 3.0;
	*/
}

static void setup_z_geometry(ham_struct* hs)
{
	const double aSI[3][3] = { { 0.28165, 0.28165, 0 }, { 0.28165, -0.28165, 0 }, { 0, 0, 0.5633 } };
	hs->L = 0.5633 / 2 / M_PI / sqrt(2); 
	setup_geometry(hs, aSI);
	/* Old code for zero strain
	 *
	const double sqrt2 = 1.4142135623730950488;

	double aa = 1 / sqrt2;
	hs->kx_vect[0] = aa;
	hs->kx_vect[1] = aa;
	hs->kx_vect[2] = 0;
	hs->ky_vect[0] = aa;
	hs->ky_vect[1] = -aa;
	hs->ky_vect[2] = 0;
	hs->kz_vect[0] = 0;
	hs->kz_vect[1] = 0;
	hs->kz_vect[2] = sqrt2;
	hs->bvect_lengths[0] = 1;
	hs->bvect_lengths[1] = 1;
	hs->bvect_lengths[2] = 1 / sqrt2;

	// Valid for zero strain only!
	hs->L = 0.5633 / 2 / M_PI / sqrt2; 
	hs->volume_factor = 1; // because the basis vectors in real space are orthogonal
	hs->volume_factor_12 = 1; // because the basis vectors in real space are orthogonal
	hs->cosine_product = 1;
	hs->cosine_product_12 = 1;
	*/
}




static void initialize_tb_geometry(ham_struct* hs, tb_geometry geometry, double exx)
{
	hs->exx = exx;
	switch (geometry)
	{
		case GEOM_XY:
			setup_xy_geometry(hs);
			break;
		case GEOM_Z:
			setup_z_geometry(hs);
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", geometry);
			abort();
			break;
	}
}

static void initialize_tb2_code(tb_geometry geometry, tb_parameterization parameterization, tb_material material)
{
	if (geometry != current_tb_geometry ||
			parameterization != current_tb_parameterization ||
			material != current_tb_material)
	{
		current_tb_geometry = geometry;
		current_tb_parameterization = parameterization;
		current_tb_material = material;
		tb2init(generate_tb2_input_filename(parameterization, geometry, material));
	}
}

/*
Common function to create (partially) bulk and plane TB hamiltonian structures.
hs - struct to be set up
dim - dimension of the hamiltonian
*/
static void setup_tightbinding_common_hamiltonian_structure(ham_struct* hs, int dim, tb_parameterization parameterization, tb_geometry geometry, tb_material material, double exx)
{
	setup_ham_struct(hs, TB_MN_X, dim, get_periodic_bc_dim(material));
	initialize_tb_geometry(hs, geometry, exx);
	hs->find_cstep = find_cstep_tb;
	hs->inverted_fermi = 1;
}

/*
 * Calculate dH/dc in exact way.
 */
void exact_tb_bulk_derivative(ham_struct* this, double c1, double c2, double c3, double dir1, double dir2, double dir3, doublecomplex matrix[][])
{
	double c[3];
	double dir[3];
	c[0] = c1;
	c[1] = c2;
	c[2] = c3;
	dir[0] = dir1;
	dir[1] = dir2;
	dir[2] = dir3;
	tb2hamilton_derivative(this->D / Delta0(this), c, dir, matrix); // D/Delta0 = magfac
}

ham_struct* create_tightbinding_bulk_hamiltonian_structure(tb_derivative derivative, tb_geometry geometry, tb_parameterization parameterization, double exx)
{
	int dim = get_bulk_dimension(geometry);

	initialize_tb2_code(geometry, parameterization, BULK);

	ham_struct* hs = alloc_memory(sizeof(ham_struct));

	setup_tightbinding_common_hamiltonian_structure(hs, dim, parameterization, geometry, BULK, exx);
	hs->p_bands_start = get_bulk_p_bands_start(parameterization, geometry);
	hs->p_bands_cnt = get_bulk_p_bands_cnt(geometry);
	hs->bril_zone_integrand = bril_zone_integrand_tb_3d;
	hs->S_z = create_S_z_bulk(dim);
	hs->S_x = create_S_x_bulk(dim);    
	hs->S_y = create_S_y_bulk(dim);    
	hs->S_minus = create_S_minus(hs->S_x, hs->S_y, dim);
	hs->S_plus = create_S_plus(hs->S_x, hs->S_y, dim);
	hs->gen_ham = calculate_tb2_bulk_hamiltonian;
	switch(derivative)
	{
		case EXACT:
			hs->gen_ham_derivative = exact_tb_bulk_derivative;
			break;
		case FIVEPOINT:
			hs->gen_ham_derivative = discretized_5pt_derivative;
			break;
		default:
			fprintf(stderr, "Unsupported derivative calculation method: %i\n", derivative);
			abort();
			break;

	}

	return hs;
}

ham_struct* create_tightbinding_plane_hamiltonian_structure(tb_derivative derivative, tb_geometry geometry, tb_parameterization parameterization, int nlayers, double Ex, double Ey, double Ez, double (*concentration_gradient)(double*), double exx)
{
	int bulk_dim = get_bulk_dimension(geometry);
	int dim = bulk_dim * nlayers;

	assert(exx == 0);

	initialize_tb2_code(geometry, parameterization, PLANE);

	plane_ham_struct* phs = alloc_memory(sizeof(plane_ham_struct));
	ham_struct* hs = (ham_struct*) phs;
	setup_tightbinding_common_hamiltonian_structure(phs, dim, parameterization, geometry, PLANE, exx);
	setup_plane_ham_struct(phs, bulk_dim, nlayers);
	phs->electric_field[0] = Ex;
	phs->electric_field[1] = Ey;
	phs->electric_field[2] = Ez;
	phs->concentration_gradient = concentration_gradient;
	hs->p_bands_start = nlayers * get_bulk_p_bands_start(parameterization, geometry);
	hs->p_bands_cnt = nlayers * get_bulk_p_bands_cnt(geometry);
	hs->bril_zone_integrand = bril_zone_integrand_tb_2d;
	double layer_thickness;
	switch (geometry)
	{
		case GEOM_XY:
			layer_thickness = 2 * M_PI; // TODO: adapt to exx != 0
			break;
		case GEOM_Z:
			layer_thickness = 2 * M_PI * sqrt(2); // TODO: adapt to exx != 0
			break;
		default:
			fprintf(stderr, "Unsupported geometry: %d\n", geometry);
			abort();
	}
	hs->width = nlayers * layer_thickness;
	hs->S_z = create_S_z_plane(bulk_dim, nlayers);
	hs->S_x = create_S_x_plane(bulk_dim, nlayers);
	hs->S_y = create_S_y_plane(bulk_dim, nlayers);
	hs->S_minus = create_S_minus(hs->S_x, hs->S_y, dim);
	hs->S_plus = create_S_plus(hs->S_x, hs->S_y, dim);
	hs->gen_ham = calculate_tb2_plane_hamiltonian;
	hs->destroy = destroy_plane_ham_struct;
	switch(derivative)
	{
		case EXACT:
			hs->gen_ham_derivative = exact_tb_plane_derivative;
			break;
		case FIVEPOINT:
			hs->gen_ham_derivative = discretized_5pt_derivative;
			break;
		default:
			fprintf(stderr, "Unsupported derivative calculation method: %i\n", derivative);
			abort();
			break;

	}

	return hs;
}

