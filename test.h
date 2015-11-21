


#pragma once
#include <afxwin.h>  // necessary for MFC to work properly
#include <math.h> 
#include <stdio.h> 
#include <queue>


#include "../../src/blepo.h"





using namespace std;
using namespace blepo;
using namespace blepo_ex;
//using namespace blepo_ex;


struct node
	{
		float x;
		float y;
		float value;
	};





int initial_filename_recognize( const CString path_input, const CString filename_input)
{
	

	if(   fopen(path_input+filename_input,"r")==0   )
	  { return 1; }
	
	else return 0;
				
}









int gaussian_kernel_compute(const float sigma_value, ImgFloat* gaussian_kernel)
{
	float f=2.5;
	int kernel_halfwidth= Round( f * sigma_value -0.5); //
	int kernel_width=2*kernel_halfwidth +1;

	//input protection
	if(kernel_width>=100 || sigma_value==0)
	{
		return 1;
		
	}


	//ImgFloat gaussian_kernel;

	(*gaussian_kernel).Reset(kernel_width,1);
	Set(gaussian_kernel,0);

	float gaussian_sum=0;
	float sum1=0;
	
	//Build Gaussian  kernel
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		(*gaussian_kernel)(ctn,0) = exp( 
						-(ctn-kernel_halfwidth) * (ctn-kernel_halfwidth)  
					   /(2 * sigma_value * sigma_value) );
		gaussian_sum+=(*gaussian_kernel)(ctn,0) ;
	
	}

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		(*gaussian_kernel)(ctn,0)/= gaussian_sum;
		sum1+=(*gaussian_kernel)(ctn,0);
	
	}


	/*printf("b.\n");
	printf("kernel_halfwidth = %d \n", kernel_halfwidth);
	printf("gaussian_sum = %f \n", gaussian_sum);
	printf("Sum1 = %f \n", sum1);
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		printf("gaussian_kernel[%d] = %f \n",ctn , (*gaussian_kernel)(ctn,0) );
	
	
	}*/



	return 0;

}



int gaussian_kernel_deriv_compute(const float sigma_value, ImgFloat* gaussian_kernel_deriv)
{
	float f=2.5;
	int kernel_halfwidth= Round( f * sigma_value -0.5); //
	int kernel_width=2*kernel_halfwidth +1;

	//input protection
	if(kernel_width>=100 || sigma_value==0)
	{
		return 1;
		
	}



	(*gaussian_kernel_deriv).Reset(kernel_width,1);
	Set(gaussian_kernel_deriv,0);

	float gaussian_deri_sum=0;
	float sum2=0;

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		(*gaussian_kernel_deriv)(ctn,0) = (ctn-kernel_halfwidth)
				*exp(  -(ctn-kernel_halfwidth) * (ctn-kernel_halfwidth)  
					   /(2 * sigma_value * sigma_value) );

		gaussian_deri_sum+=(*gaussian_kernel_deriv)(ctn,0) * ctn;
	
	}

	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		(*gaussian_kernel_deriv)(ctn,0)= (-1)*(*gaussian_kernel_deriv)(ctn,0)/gaussian_deri_sum;
		sum2+=(*gaussian_kernel_deriv)(ctn,0);
	
	}

	
	/*printf("\n\n");
	printf("gaussian_deri_sum = %f \n", gaussian_deri_sum);
	printf("Sum2 = %f \n", sum2);
	printf("With flip \n");
	for(int ctn=0; ctn<kernel_width ; ctn++)
	{
		
		printf("gaussian_deri_kernel[%d] = %f \n",ctn , (*gaussian_kernel_deriv)(ctn,0));
	}*/

	return 0;
}




void separable_convolution_x( const ImgFloat &img_in,  const int kernel_width, const ImgFloat gaussian_kernel, ImgFloat* img_out ) 
{


	(*img_out).Reset( img_in.Width(),img_in.Height() );
	Set(img_out, 0);

	
	for( int y=0; y< img_in.Height(); y++ ) 
	{
		for( int x=(kernel_width-1)/2; x< img_in.Width()-(kernel_width-1)/2; x++ ) 
		{
			float tmp_value=0;
			for(int ctn=0; ctn<kernel_width ; ctn++)
			{

				if(	  ( (x+ (kernel_width-1)/2 - ctn) < 0 ) 
					||( (x+ (kernel_width-1)/2 - ctn) >= img_in.Width() )   
				  )

				{
					tmp_value+=gaussian_kernel(ctn,0) * 0.0f;				
				}

				else 
				{
					tmp_value+=gaussian_kernel(ctn,0) * img_in(x+ (kernel_width-1)/2 - ctn,y);
				}

			}//for ctn
			//if( tmp_value>255 ) tmp_value=255;
			if( tmp_value<0 ) tmp_value=0;
			(*img_out)(x,y)=tmp_value;
			

		}//for x
	}//for y

	//return 0;

}



//convolution_y   with kernel
void separable_convolution_y( const ImgFloat &img_in,  const int kernel_width, const ImgFloat gaussian_kernel, ImgFloat* img_out ) 
{


	(*img_out).Reset( img_in.Width() , img_in.Height() );
	Set(img_out, 0);

	for( int x=0; x< img_in.Width(); x++) 
	{
		for(  int y=(kernel_width-1)/2; y< img_in.Height()-(kernel_width-1)/2; y++ ) 
		{
			float tmp_value=0;

			for(int ctn=0; ctn<kernel_width ; ctn++)
			{

				if(	  ( (y+ (kernel_width-1)/2 - ctn) < 0 ) 
					||( (y+ (kernel_width-1)/2 - ctn) >= img_in.Height() )   
				  )

				{
					tmp_value+=gaussian_kernel(ctn,0) * 0.0f;				
				}

				else 
				{
					tmp_value+=gaussian_kernel(ctn,0) * img_in(x,y + (kernel_width-1)/2 - ctn );
				}

			}//for ctn
			
			//if( tmp_value>255 ) tmp_value=255;
			if( tmp_value<0 ) tmp_value=0;
			(*img_out)(x,y)=tmp_value;
			

		}//for x
	}//for y

	//return 0;

}



//convolution x and y  with kernel
int separable_convolution_2d( const ImgGray &img_in_g,  const int kernel_width,const ImgFloat gaussian_kernel ,const ImgFloat gaussian_deri_kernel, ImgFloat* img_out, const int mode=0 ) 
{
	ImgFloat img_in_f;
	ImgFloat img_mid;
	

	Convert(img_in_g, &img_in_f);


	if(mode == 0) // x  y'
	{
		separable_convolution_y( img_in_f,  gaussian_deri_kernel.Width(), gaussian_deri_kernel, &img_mid );
		separable_convolution_x( img_mid,  gaussian_kernel.Width(), gaussian_kernel,  img_out);
		
		return 0;
		
	}

	else if(mode == 1)// y x'
	{
		separable_convolution_x( img_in_f,  gaussian_deri_kernel.Width(), gaussian_deri_kernel, &img_mid );
		separable_convolution_y( img_mid,  gaussian_kernel.Width(), gaussian_kernel,  img_out);
		
		return 0;
	}

	else
	{
		return 1;
	}
		
}

void Gradient_gaussian(const float sigma_grad, const ImgGray &img_in_g
					   , ImgFloat* img_grad_y, ImgFloat* img_grad_x, 
					   ImgFloat* img_grad_xx, ImgFloat* img_grad_yy, ImgFloat* img_grad_xy
					   ,int* half_Width )
{
	
	ImgFloat gaussian_kernel;
	ImgFloat gaussian_kernel_deriv;
	
	
	if(gaussian_kernel_compute(sigma_grad, &gaussian_kernel) ==1) printf("Error:gk\n"); 
	if(gaussian_kernel_deriv_compute(sigma_grad, &gaussian_kernel_deriv) ==1) printf("Error:gdk\n"); 
	*half_Width=gaussian_kernel.Width()/2;
	//Grad y
	if( separable_convolution_2d( img_in_g,  gaussian_kernel.Width(), gaussian_kernel , gaussian_kernel_deriv, img_grad_y, 0 )  ==1) printf("gy\n"); 
	
	//Grad x
	if( separable_convolution_2d( img_in_g,  gaussian_kernel.Width(), gaussian_kernel , gaussian_kernel_deriv, img_grad_x, 1 )   ==1) printf("gx\n"); 

	ImgFloat img_in_f; Convert(img_in_g, &img_in_f);
	Gradient(img_in_f, sigma_grad, img_grad_x, img_grad_y);

	Multiply(*img_grad_x, *img_grad_x, img_grad_xx);
	Multiply(*img_grad_y, *img_grad_y, img_grad_yy);
	Multiply(*img_grad_x, *img_grad_y, img_grad_xy);


}





void Compute_window_sum(const int window_size, ImgFloat* img_io_grad_xx, ImgFloat* img_io_grad_yy, ImgFloat* img_io_grad_xy)
{
  ImgFloat img_mid;
  ImgFloat kernel_x(window_size,1), kernel_y(1,window_size);
  Set(&kernel_x, 1);
  Set(&kernel_y, 1);
  Convolve(*img_io_grad_xx, kernel_x, &img_mid);
  Convolve(img_mid, kernel_y, img_io_grad_xx);
  Convolve(*img_io_grad_yy, kernel_x, &img_mid);
  Convolve(img_mid, kernel_y, img_io_grad_yy);
  Convolve(*img_io_grad_xy, kernel_x, &img_mid);
  Convolve(img_mid, kernel_y, img_io_grad_xy);
}






void capture_features(const ImgGray& img_in_g, queue<node>* points_featured_sec )
{
	//while (!(*points_featured_sec).empty()) (*points_featured_sec).pop();
	ImgFloat img_grad_x;
	ImgFloat img_grad_y;

	ImgFloat img_grad_xx;
	ImgFloat img_grad_yy;
	ImgFloat img_grad_xy;

	int grad_half_Width=0;

	//test sigma
	const float value_sigma_grad=1.5f;
	Gradient_gaussian(value_sigma_grad, img_in_g
					   , &img_grad_y, &img_grad_x, 
					   &img_grad_xx, &img_grad_yy, &img_grad_xy 
					   ,&grad_half_Width);
	
	int window_size=3;
    Compute_window_sum(window_size, &img_grad_xx, &img_grad_yy, &img_grad_xy);


	ImgFloat img_min_eigvalue_map;
	

	//data structure
	node points_featured;
	queue<node> points_sq_featured;
	int points_sq_deepth=0;

	(img_min_eigvalue_map).Reset(img_grad_x.Width(), img_grad_x.Height());
	Set(&img_min_eigvalue_map, 0);
	

	
	
	for (int y=grad_half_Width ; y<img_grad_x.Height()-grad_half_Width ; y++)
	{
		for (int x=grad_half_Width ; x<img_grad_x.Width()-grad_half_Width ; x++)
		{
			float pixel_value_g_xx = img_grad_xx(x, y);
			float pixel_value_g_yy = img_grad_yy(x, y);
			float pixel_value_g_xy = img_grad_xy(x, y);
			
			float min_eig_for_pixel;
			min_eig_for_pixel = (pixel_value_g_xx + pixel_value_g_yy - sqrt((pixel_value_g_xx - pixel_value_g_yy)*(pixel_value_g_xx - pixel_value_g_yy) + 4.0f*pixel_value_g_xy*pixel_value_g_xy))/2.0f;
			
			(img_min_eigvalue_map)(x, y) = min_eig_for_pixel;
			
			/*points_featured.x=(float)x;
			points_featured.y=(float)y;
			points_featured.value=min_eig_for_pixel;
			points_sq_featured.push(points_featured);*/
			//points_sq_deepth++;

		}
	}


	//threshold
	const float min_eig = 75.0f;
	ImgFloat img_threshold_min_eig_map;
	(img_threshold_min_eig_map).Reset(img_grad_x.Width(), img_grad_x.Height());
	Set(&img_threshold_min_eig_map, 0);

	float *pixel_min_eig_map= (img_min_eigvalue_map).Begin();
	float *pixel_th_min_eig_map= (img_threshold_min_eig_map).Begin();


	while(pixel_min_eig_map!=(img_min_eigvalue_map).End() )
	{
		if(	*pixel_min_eig_map> min_eig )
		{
			*pixel_th_min_eig_map=*pixel_min_eig_map;		
		}
	
		pixel_min_eig_map++;
		pixel_th_min_eig_map++;	
	}


	//NMS
	ImgFloat img_nms_min_eig_map;
	(img_nms_min_eig_map).Reset(img_grad_x.Width(), img_grad_x.Height());
	Set(&img_nms_min_eig_map, 0);
	//img_nms_min_eig_map=img_threshold_min_eig_map;

	for (int y=grad_half_Width ; y<img_grad_x.Height()-grad_half_Width ; y++)
	{
		for (int x=grad_half_Width ; x<img_grad_x.Width()-grad_half_Width ; x++)
		{
			float comp_value[5]={0};
			comp_value[1] = (img_threshold_min_eig_map)(x-1,y)>(img_threshold_min_eig_map)(x+1,y)? (img_threshold_min_eig_map)(x-1,y): (img_threshold_min_eig_map)(x+1,y);
			comp_value[2] = (img_threshold_min_eig_map)(x-1,y-1)>(img_threshold_min_eig_map)(x+1,y+1)?  (img_threshold_min_eig_map)(x-1,y-1):(img_threshold_min_eig_map)(x+1,y+1);
			comp_value[3] =(img_threshold_min_eig_map)(x,y-1)>(img_threshold_min_eig_map)(x,y+1)? (img_threshold_min_eig_map)(x,y-1):(img_threshold_min_eig_map)(x,y+1);
			comp_value[4] =(img_threshold_min_eig_map)(x+1,y-1)> (img_threshold_min_eig_map)(x-1,y+1)?(img_threshold_min_eig_map)(x+1,y-1): (img_threshold_min_eig_map)(x-1,y+1);
			
			bool sign_max=false;
			for(int ctn=1;ctn<5;ctn++)
			{	
				if(  (img_threshold_min_eig_map)(x,y)>comp_value[ctn] ) sign_max=true;
				else sign_max=false;						
			}
			if(sign_max==true) 
			{
					(img_nms_min_eig_map)(x,y)=(img_threshold_min_eig_map)(x,y);
					points_featured.x=(float)x;
					points_featured.y=(float)y;
					points_featured.value=(img_threshold_min_eig_map)(x,y);
					points_sq_featured.push(points_featured);

			}
		
		}
	}


	int filled_window=3;
	ImgGray img_featured(img_grad_x.Width(), img_grad_x.Height());
	Set(&img_featured, 0);


	//queue<node> points_featured_sec;

	while( !(points_sq_featured).empty())
	{
		const node current_point = points_sq_featured.front();
		points_sq_featured.pop();
		if (current_point.value < min_eig)  continue;
		
		if (!img_featured((int)current_point.x , (int)current_point.y)  )  
		{
			(*points_featured_sec).push(current_point);
			int left   = Max((int)current_point.x - filled_window, 3);
			int right  = Min((int)current_point.x + filled_window, img_grad_x.Width()-3);
		
			int up    = Max((int)current_point.y - filled_window, 3);
			int down  = Min((int)current_point.y + filled_window, img_grad_x.Height()-3);
			Set(&img_featured, Rect(left, up, right, down), 1);
			//printf("yy ");
		}
	}

	/*ImgBgr img_wt_featured;
	Convert(img_in_g, &img_wt_featured);
	int queue_size=(*points_featured_sec).size();
	for(int ctn=0;ctn<=queue_size;ctn++)
	{
	
		node point = (*points_featured_sec).front();
		(*points_featured_sec).pop();
		DrawDot(Point(point.x, point.y), &img_wt_featured, Bgr(0,0,255));
		(*points_featured_sec).push(point);
		printf("xx ");
	}
	Figure fig_wt_featured;
	fig_wt_featured.SetTitle("image_with_featured");
	fig_wt_featured.Draw(img_wt_featured);*/


	/*Figure fig_min_eigvalue_map;
	fig_min_eigvalue_map.SetTitle("img_min_eigvalue_map");
	fig_min_eigvalue_map.Draw(*img_min_eigvalue_map);
 
	Figure fig_threshold_min_eig_map;
	fig_threshold_min_eig_map.SetTitle("img_threshold_min_eig_map");
	fig_threshold_min_eig_map.Draw(*img_threshold_min_eig_map);*/

	/*Figure fig_nms_min_eig_map;
	fig_nms_min_eig_map.SetTitle("img_nms_min_eig_map");
	fig_nms_min_eig_map.Draw(img_nms_min_eig_map);
*/
	


}




int tracking_pair(const ImgGray& img_in_1, const ImgGray& img_in_2
				   ,queue<node>* points_featured 
				   ,const float sigma_value
				   ,const int window_size)
				   //,const ImgFloat gaussian_kernel, const ImgFloat gaussian_kernel_deriv)
{

	//release queue2
	//while (!(*points_featured_2).empty()) (*points_featured_2).pop();
	if( (*points_featured).empty()==true ) return 1; 


	ImgFloat img_in_1_f;
	ImgFloat img_in_2_f;
	Convert(img_in_1,&img_in_1_f);
	Convert(img_in_2,&img_in_2_f);

	ImgFloat img_grad_x;
	ImgFloat img_grad_y;

	ImgFloat img_grad_xx;
	ImgFloat img_grad_yy;
	ImgFloat img_grad_xy;

	int kernel_halfwidth = 0;

	Gradient_gaussian(sigma_value, img_in_1
					   , &img_grad_y, &img_grad_x
					   , &img_grad_xx, &img_grad_yy, &img_grad_xy
					   ,&kernel_halfwidth );
	
	//int window_size=3;
    Compute_window_sum(window_size, &img_grad_xx, &img_grad_yy, &img_grad_xy);



	const int window_size_half= window_size/2; 
	int queue_size_pt_featured=(*points_featured).size();
	for (int ctn=0 ; ctn<queue_size_pt_featured; ctn++)
	{
		node current_point = (*points_featured).front();
		(*points_featured).pop();

		if( (int)current_point.x-kernel_halfwidth<window_size_half 
			|| (int)current_point.y-kernel_halfwidth<window_size_half 
			|| (int)current_point.x+window_size_half+kernel_halfwidth+1>=img_in_1.Width() 
			||(int)current_point.y+window_size_half+kernel_halfwidth+1>=img_in_1.Height() 
		  ) continue;

		/*float du_total = 0;
		float dv_total = 0;*/

		//Z matrix
		float g_xx = img_grad_xx(Round(current_point.x) , Round(current_point.y) );
		float g_yy = img_grad_yy(Round(current_point.x) , Round(current_point.y) );
		float g_xy = img_grad_xy(Round(current_point.x) , Round(current_point.y) );

		// windows in two img
		Rect region(Round(current_point.x-window_size_half), Round(current_point.y-window_size_half)
					,Round(current_point.x+window_size_half+1) , Round(current_point.y+window_size_half+1) );
		
		ImgFloat img_in_1_region;
		Extract(img_in_1_f, region, &img_in_1_region);

		ImgFloat img_in_2_region;
		Extract(img_in_2_f, region, &img_in_2_region);

		ImgFloat img_grad_x_region;
		Extract(img_grad_x, region, &img_grad_x_region);
	
		ImgFloat img_grad_y_region;
		Extract(img_grad_y, region, &img_grad_y_region);

		// error vector
		float error_vector_x=0, error_vector_y=0;

		// du dv
		float du_sum=0;
		float dv_sum=0;
		
		//u

		//exit point
		int iter=0;
		int max_iter = 10;
		float min_offset =0.1f;

		bool loop=false;
		//int u=0;
		while (!loop)
		{
			//find e
			if(IsSameSize(img_in_1_region, img_in_2_region) 
				&& IsSameSize(img_in_1_region, img_grad_x_region) 
				&& IsSameSize(img_in_1_region, img_grad_y_region)  )
			{
				const float* point_1 = img_in_1_region.Begin();
				const float* point_2 = img_in_2_region.Begin();
				const float* point_x = img_grad_x_region.Begin();
				const float* point_y = img_grad_y_region.Begin();

				error_vector_x = 0;
				error_vector_y = 0;

				while (point_1 != img_in_1_region.End())
				{

					error_vector_x += (*point_x)* (*point_1 - *point_2);
					error_vector_y += (*point_y)* (*point_1 - *point_2);

					point_1++;
					point_2++;
					point_x++;
					point_y++;
				}
			}



			// solve Zd = e
			float determinant = g_xx * g_yy - g_xy * g_xy; 
			float du = (g_yy * error_vector_x - g_xy * error_vector_y) / determinant;
			float dv = (g_xx * error_vector_y - g_xy * error_vector_x) / determinant;
			du_sum += du;
			dv_sum += dv;

			// shift image and recompute 'e'
			
				
				img_in_2_region.Reset(window_size, window_size);
				for (int y = 0 ; y < window_size; y++)
				{
					for (int x = 0 ; x < window_size ; x++)
					{
						float xf=x+current_point.x+du_sum-window_size_half;
						float yf=y+current_point.y+dv_sum-window_size_half;

						if (xf<0)  xf = 0.0f;
						else if (xf >= img_in_2.Width()-1)  xf = img_in_2.Width() - 1.01f;
						if (yf<0)  yf = 0.0f;
						else if (yf >= img_in_2.Height()-1)  yf = img_in_2.Height() - 1.01f;
						
						int x0 = (int)xf;
						int y0 = (int)yf;

						float xa = xf - x0;
						float yb = yf - y0;
												
						img_in_2_region(x, y) =  (1-xa)*(1-yb)*img_in_2_f(x0, y0) + xa*(1-yb)*img_in_2_f(x0+1, y0)
									+ (1-xa)*yb*img_in_2_f(x0,y0+1) + xa*yb*img_in_2_f(x0+1, y0+1) ;
					}
				}
			

			//for loop 
			iter++;
			if(iter > max_iter || ((fabs(du) < min_offset) && (fabs(dv) < min_offset)) ) loop=true;
		}
		current_point.x += du_sum;
		current_point.y += dv_sum;
		(*points_featured).push(current_point);
			
	}
	
	/*ImgBgr img_wt_featured;
	Convert(img_in_2, &img_wt_featured);
	int queue_size_pt_f=(*points_featured).size();
	for(int ctn=0;ctn<=queue_size_pt_f;ctn++)
	{
	
		node point = (*points_featured).front();
		(*points_featured).pop();
		DrawDot(Point(point.x, point.y), &img_wt_featured, Bgr(0,0,255));
		(*points_featured).push(point);
		printf("xx ");
	}
	Figure fig_wt_featured;
	fig_wt_featured.SetTitle("image_with_featured");
	fig_wt_featured.Draw(img_wt_featured);*/

	return 0;
}





//int tracking_sequence( const ImgGray img_loaded_sequence_current_g[], queue<node>* points_featured[]
//				   
//				   ,const int max_sequence_num,const int window_size
//				   ,const ImgFloat gaussian_kernel, const ImgFloat gaussian_kernel_deriv )
//{
//
//	for(int ctn=0; ctn<max_sequence_num;ctn++)
//	{
//		if (tracking_pair(img_loaded_sequence_current_g[ctn], img_loaded_sequence_current_g[ctn+1]
//				   ,points_featured[ctn]
//				   ,window_size
//				   ,gaussian_kernel, gaussian_kernel_deriv)==0 )   printf("get");
//
//	}
//
//	return 0;
//}










