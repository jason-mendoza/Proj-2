#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color= new pixel[width * height];
    state.image_depth = new float[width * height];
    int length = height * width;
    for(int i = 0; i < length; i++){
	state.image_color[i] = make_pixel(0,0,0); 
	state.image_depth[i] = 2;
    }
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{

	for(int i = 0; i < state.num_vertices; i +=3){
		data_geometry ** three = new data_geometry*[3];
		for(int q = 0; q < 3; q++){
			three[q] = new data_geometry;
		}
		for(int q  = 0; q < 3; q++){
			three[q]->data = new float[MAX_FLOATS_PER_VERTEX];
		}
		
		for(int j = 0; j < state.floats_per_vertex; j++){
			three[0]->data[j] = state.vertex_data[j + state.floats_per_vertex * (i) ];
			three[1]->data[j] = state.vertex_data[j + state.floats_per_vertex * (i + 1) ];
			three[2]->data[j] = state.vertex_data[j + state.floats_per_vertex * (i + 2) ];
		}
		
		//data_geometry** out = new data_geometry*[3];
		for(int j = 0; j < 3; j++){
			data_vertex temp;
			temp.data = three[j]->data;
			state.vertex_shader((const data_vertex)temp,*three[j], state.uniform_data);
		}
		//std::cout << three[0] -> data[0] << std::endl;
		//for(int p = 0; p < 3; p++){
		//	three[p]->gl_Position /= three[p]->gl_Position[3];
		//}
		rasterize_triangle(state, (const data_geometry **) three);
		delete[] three[0] -> data; delete[] three[1] -> data; delete[] three[2]->data;
		delete three[0]; delete three[1]; delete three[2]; delete[] three;
	}
/*
    int q  = 0;
    data_geometry * in[3] = {0,0,0};
    float * data_temp = 0;
    for( int i = 0; i < (state.num_vertices * state.floats_per_vertex); i++){
	if(  state.floats_per_vertex % i  == 0 ){
		data_temp = new float[MAX_FLOATS_PER_VERTEX];
		in[q] = new data_geometry;
	}
		//float * data2_temp = new float[MAX_FLOATS_PER_VERTEX];
		//float * data3_temp = new float[MAX_FLOATS_PER_VERTEX];
	data_temp[q] = state.vertex_data[i];
	if(state.floats_per_vertex % (i + 1) == 0 && i != 0){
		in[q]->data = data_temp;
		delete data_temp;
		q++;
		if( (state.floats_per_vertex * 3) % (i + 1) == 0){ 
			const data_geometry * place[3] = {in[0], in[1], in[2]};
			rasterize_triangle(state, place);
			delete in[0]; delete in[1]; delete in[2];
			in[0] = 0; in[1] = 0; in[2] = 0;
			q = 0;
		}	
	}	
	
    }
*/
    
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    int i;
    int j;
    int width = state.image_width;
    int height = state.image_height;

    int ax;
    int ay;
    int bx;
    int by;
    int cx;
    int cy;

    float AreaABC;
    float AreaPBC;
    float AreaAPC;
    float AreaABP;
    float alpha;
    float beta;
    float gamma;
    unsigned int image_index = 0;
    //data_geometry* out = (data_geometry*)in;
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    /*
    for(int index = 0; index < 3; index++){
	//data_geometry* out = new data_geometry[3];
    	//state.vertex_shader(temp, out[index], state.uniform_data);
		//for(int q = 0; q < 2; q++){
		//	out[index].gl_Position[q] /= out[index].gl_Position[3];
		//}
	i = width / 2.0 * in[index]->gl_Position[0] + width/2.0 - (0.5); 
	j = height / 2.0 * in[index]->gl_Position[1] + height /2.0 - (0.5);

	image_index = i + j * width;
	state.image_color[image_index] = make_pixel(255,255,255);
    }
    */
    ax = (width / 2.0) * (in[0]->gl_Position[0] / in[0]->gl_Position[3]) + (width/2.0) - (0.5);
    ay = (height / 2.0) * (in[0]->gl_Position[1] / in[0]->gl_Position[3]) + (height/ 2.0) - (0.5) ;
    bx = (width / 2.0) * (in[1]->gl_Position[0] / in[1]->gl_Position[3]) + (width/2.0) - (0.5);
    by = (height / 2.0) *(in[1]->gl_Position[1] / in[1]->gl_Position[3]) + (height / 2.0) - (0.5);
    cx = (width / 2.0) * (in[2]->gl_Position[0] / in[2]->gl_Position[3]) + (width / 2.0) - (0.5);
    cy = (height / 2.0) * (in[2]->gl_Position[1] / in[2]->gl_Position[3]) + (height / 2.0) - (0.5);
    
    AreaABC = 0.5 * (((bx * cy) - (cx * by)) - ((ax * cy)-(cx * ay)) + ((ax * by)-(bx * ay)));   
    for(int px = 0; px < width; px++){
    	for(int py = 0; py < height; py++){
		AreaPBC = 0.5 * (((bx * cy) - (cx * by)) + ((by-cy)*px) + ((cx -bx)*py));
		AreaAPC = 0.5 * (((cx * ay) - (ax * cy)) + ((cy - ay) * px) + ((ax - cx) * py) );  
		AreaABP = 0.5 * (((ax * by) - (bx*ay)) + ((ay-by)*px) + ((bx-ax)*py));
		
		alpha = AreaPBC /  AreaABC;
		beta = AreaAPC / AreaABC;
		gamma = AreaABP / AreaABC;
		
		image_index =  px + py * width;
		if(alpha >= 0 && beta >= 0 && gamma >= 0){
			//state.image_color[image_index] = make_pixel(255,255,255);
			float alpha_p;
			float beta_p;
			float gamma_p;
			//int index = px + py * state.image_width;
			auto *data = new float[MAX_FLOATS_PER_VERTEX];
			data_fragment f{data};
			data_output frag_out;

			float depth1 = alpha * in[0]->gl_Position[2] + beta * in[1]->gl_Position[2] + gamma * in[2]->gl_Position[2];
			if(depth1 > state.image_depth[px + py * width]){
				continue;
			}

			for(int q = 0; q < state.floats_per_vertex; q++){
				float k_gour;
				switch(state.interp_rules[q]){
					case interp_type::flat:
					f.data[q] = in[0]->data[q];
					break;
					
					case interp_type::smooth:
					k_gour = (alpha/in[0] -> gl_Position[3] + 
							beta/in[1] ->gl_Position[3] + 
							gamma/in[2]->gl_Position[3]);
					alpha_p = alpha /  (k_gour * (in[0] -> gl_Position[3]) );
					beta_p = beta /  (k_gour * (in[1] -> gl_Position[3]));
					gamma_p = gamma / (k_gour * (in[2] -> gl_Position[3]) );
					f.data[q] = (alpha_p * in[0]->data[q]) + (beta_p * in[1]->data[q]) + (gamma_p * in[2] ->data[q]);
					break;
					
					case interp_type::noperspective:
					f.data[q] = alpha*in[0]->data[q] + beta*in[1]->data[q] + gamma*in[2]->data[q];
					break;
					
					default:
					break;						
				}
			}
			state.fragment_shader( f, frag_out, state.uniform_data );
			state.image_color[px + py * width] = make_pixel(frag_out.output_color[0] * 255,
								 frag_out.output_color[1] * 255,
								 frag_out.output_color[2] * 255);
			state.image_depth[px + py * width] = depth1;
		} 
	}
    }

}
