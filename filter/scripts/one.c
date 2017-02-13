/*
Fast Artificial Neural Network Library (fann)
Copyright (C) 2003 Steffen Nissen (lukesky@diku.dk)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include "floatfann.h"

int main(int argc, char ** argv)
{
	fann_type *calc_out;
	fann_type input[2];

	struct fann *ann = fann_create_from_file(argv[2]);
        struct fann_train_data *data = fann_read_train_from_file(argv[1]);
	
	unsigned int i = 0;
	int j = 0;
	unsigned int decimal_point;
	for(i = 0; i < data->num_data; i++){
		calc_out = fann_run(ann, data->input[i]);
		printf("%i\n",i+1);
		for(j=0;j<1;j++)
		{
		  printf(" %i %f  %f\n",
			 j+1,
			 data->output[i][j],
			 calc_out[j]);
			 
		}

	}
	
	fann_destroy(ann);
	return 0;
}
