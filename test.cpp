



//command parameters: statue_seq/img0%d.bmp  588 618    0.5  11

#include <afxwin.h>  // necessary for MFC to work properly
//#include <afxpriv.h> //for W2A
#include <locale.h>  
//#include <stdio.h>  

#include <math.h> 
#include "test.h"
#include "../../src/blepo.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif





using namespace blepo;
using namespace blepo_ex;


int main(int argc, const char* argv[], const char* envp[])
{
	
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", hModule);
		return 1;
	}
	

	printf("argc = %d\n", argc);
	for (int i=0; i<argc ; i++) 
	{
		printf("argv[%d]=%s\n\n", i, argv[i]);
	}
	
	if( argc < 5)	
	{
		printf("\nError: lack enough command parameters \n");
		printf("Act: Program halt, Please close \n");
		EventLoop();
		return 0;
				
	}
	else if(argc > 6)
	{
		printf("\nWarning: more command parameters than use \n");
	}
	

	/*
	a.	Reads 5 command-line parameters,   
	*/
	printf(_T("\nNote: Example of Commandline Parameters:\n  statue_seq/img0%d.bmp  588 618    0.5  11\n"),"%d" );
	CString path_blepo_images ="../../images/";
	//CString path_blepo_images_sq ="../../images/flowergarden/";



	CString filename_format=_T(argv[1]);//loading file_format
	CString filename_first_frame=argv[2];//loading first frame
	CString filename_last_frame=argv[3];//loading last_frame
	CString filename_value_sigma=argv[4];//loading sigma
	CString filename_value_window_size=argv[5];//loading window size

	printf("Note: Input filename_format is %s \n",filename_format);
	int value_first_frame=_ttoi(filename_first_frame);
	printf("Note: Input value_first_frame is %d \n",value_first_frame);

	int value_last_frame=_ttoi(filename_last_frame);
	printf("Note: Input value_last_frame is %d \n",value_last_frame);

	 

	float value_sigma=(float) _ttof(filename_value_sigma);
	printf("Note: Input value_sigma is %f \n",value_sigma);

	int value_window_size=_ttoi(filename_value_window_size);
	printf("Note: Input value_window_size is %d \n",value_window_size);


	int value_sequence_current=value_first_frame;

	//defence
	
	if(value_last_frame-value_first_frame<0 || value_last_frame-value_first_frame>99)
	{	
		printf("Error: too large or wrong sequence in argv[2],argv[3]");
		EventLoop();
		return 0;
	}

	if(value_sigma<0 || value_sigma>10)
	{	
		printf("Error: unsafe value in argv[4]");
		EventLoop();
		return 0;
	}

	if(value_window_size<0 || value_window_size>50)
	{	
		printf("Error: unsafe value in argv[5]");
		EventLoop();
		return 0;
	}
	

	int max_sequence_num=value_last_frame-value_first_frame;
	ImgGray img_loaded_sequence_current_g[100];
	//ImgFloat img_nms_min_eig_map[100];
	queue<node> points_featured[100];



	for( int ctn= value_first_frame; ctn<=value_last_frame; ctn++)
	{
		CString filename_sequence_current;
		filename_sequence_current.Format(_T(filename_format), ctn);// ‘filename’ is now “img0%d.bmp”
		CString filename_forloading_current=path_blepo_images+filename_sequence_current;

		int status_filename_loading_error=0;
		int status_filename_current=initial_filename_recognize(path_blepo_images, filename_sequence_current);
		//printf("%s%s \n",path_blepo_images ,filename_sequence_current );


		if( status_filename_current==1 )
		{
			status_filename_loading_error=1;
			printf("\nError: File cannot be found : %s%s \n",path_blepo_images ,filename_sequence_current );
			EventLoop();
			return 0;
			//break;	
		}
		
	 

		//Figure and title(put out)
		/*Figure fig_loaded_current;
		fig_loaded_current.SetTitle("img_loaded_sequence_current");*/
		//load and convert left filename
		Load(filename_forloading_current, &img_loaded_sequence_current_g[ctn-value_first_frame]);
		
	}
	capture_features(img_loaded_sequence_current_g[0], &points_featured[0] );


	printf("C done\n"); 
	if( (points_featured[0]).empty()==true ) 	printf("C done1\n"); 
 

	ImgBgr img_wt_featured;
	Convert(img_loaded_sequence_current_g[0], &img_wt_featured);
	int queue_size_pt_f=points_featured[0].size();
	for(int ctn=0;ctn<queue_size_pt_f;ctn++)
	{
	
		node point = points_featured[0].front();
		points_featured[0].pop();
		DrawDot(Point((int)point.x, (int)point.y), &img_wt_featured, Bgr(0,0,255));
		points_featured[0].push(point);
		//printf("xx ");
	}
	Figure fig_wt_featured;
	fig_wt_featured.SetTitle("First_image_with_featured");
	fig_wt_featured.Draw(img_wt_featured);

	
	

	ImgBgr img_wt_featured_motion;
	Figure fig_wt_featured_motion;
	fig_wt_featured_motion.SetTitle("image_with_featured_motion");


	for(int ctn=0; ctn<max_sequence_num;ctn++)
	{
		if (tracking_pair(img_loaded_sequence_current_g[ctn], img_loaded_sequence_current_g[ctn+1]
				   ,&points_featured[0]
				   ,value_sigma
				   ,value_window_size)==1 )   
		//TrackFeatures(&points_featured[0], img_loaded_sequence_current_g[ctn], img_loaded_sequence_current_g[ctn+1]);
		printf("cpp%d\n",ctn);


		
		Convert(img_loaded_sequence_current_g[ctn+1], &img_wt_featured_motion);
		int queue_size_pt_f=(points_featured)[0].size();
		for(int i=0;i<=queue_size_pt_f;i++)
		{

			node point = (points_featured)[0].front();
			(points_featured)[0].pop();
			DrawDot(Point((int)point.x, (int)point.y), &img_wt_featured_motion, Bgr(0,0,255));
			(points_featured)[0].push(point);
			//printf("xx ");
		}
		fig_wt_featured_motion.Draw(img_wt_featured_motion);


		Sleep(300); 
		
	}


	printf("D done\n"); 
	

	EventLoop();
	return 0;
}










