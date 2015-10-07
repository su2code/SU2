/*
 *     Copyright (C) 2012  Adam Jirasek
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     contact: libm3l@gmail.com
 * 
 */



/*
 *     Function Client1.c
 *
 *     Date: 2013-02-23
 * 
 * 
 *     Description:
 * 
 *
 *     Input parameters:
 * 
 *
 *     Return value:
 * 
 * 
 *
 *     Modifications:
 *     Date		Version		Patch number		CLA 
 *
 *
 *     Description
 * 
 */


#include "libm3l.h"
#include "lsipdx.h"

int main(int argc, char *argv[])
{
	node_t *Gnode=NULL, *Snode=NULL, *FoundNode=NULL, *TmpNode=NULL;
	size_t i, niter, dim[1], tot_dim, narray;

	lmint_t sockfd, portno;

        socklen_t clilen;
        struct sockaddr_in cli_addr;
	lmchar_t *name ="CFD2SIM";
	lmchar_t *name1="SIM2CFD";

	lmdouble_t *P, dy, *tmpfloat, *x, *y, *z, *time, sign;
	lmdouble_t psi, theta, phi;
	
	find_t *SFounds;
	
	opts_t opts, *Popts_1;
	
	FILE *fp;
	
	
	client_fce_struct_t InpPar, *PInpPar;

	PInpPar = &InpPar;
/*
 * get port number
 */
	if (argc < 3) {
		fprintf(stderr,"ERROR, no IPaddress and port number provided\n");
		exit(1);
	}
 	portno = atoi(argv[2]);
/*
 * open socket - because we use more then just send - receive scenario
 * we need to open socket manualy and used Send_receive function with hostname = NULL, ie. as server
 * portno is then replaced by socket number
 */
	niter = 0;
 	while(1){

 		printf("\n\n--------------------------------    i = %ld\n\n", ++niter);
/*
 * open socket
 */
		PInpPar->channel_name = name;
		PInpPar->SR_MODE = 'R';
		if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
			Error("wrong client mode");
		Popts_1 = &opts;
		m3l_set_Send_receive_tcpipsocket(&Popts_1);
	
		if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
			Error("client_sender: Error when opening socket");
		
		
		Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);
		
		printf(" Data from Edge received\n");
	
		if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
			Error("CatData");
// 		if(m3l_Cat(Gnode, "--detailed", "-P", "-L",  "*",   (char *)NULL) != 0)
// 		Error("CatData");

	if( (SFounds = m3l_Locate(Gnode, "/CFD_2_SIM/ForcesMoments", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_edge2out: More then one DX data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_edge2out: Did not find 1st data pointer");
		tot_dim = m3l_get_List_totdim(FoundNode);
		narray = tot_dim;
		
		printf(" Size of array is %ld\n", narray);
		
		if( (x = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_edge2out: Did not find DX data pointer");
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_edge2out: P not found\n");
	}
	
	SFounds = m3l_Locate(Gnode, "/CFD_2_SIM/Time", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	time = (lmdouble_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);

	printf("Time is %lf  tot-dim is %ld\n", *time, tot_dim);
	
	
/*
 *  close socket
 */
	if( close(sockfd) == -1)
		Perror("close");
	
// 	printf(" ENTER \n");
// 	getchar();
	
	
// 	sleep(10);

	
	
/*
 * open socket for sending data back
 */
	
	
	PInpPar->channel_name = name1;
	PInpPar->SR_MODE = 'S';
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("wrong client mode");
	Popts_1 = &opts;
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	
	if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");

	if(  (Snode = m3l_Mklist("SIM_2_CFD", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("m3l_Mklist");
		
// 	dim[0] = 9;
	dim[0] = 3;
/*
 * add iteration number
 */
	if(  (TmpNode = m3l_Mklist("Angles", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
/*
 * 
 *  psi - yaw
 *  theta - pitch
 *  phi - roll
 */
	psi = *time*2;
	theta=*time*4;
	phi = *time*22.5;
/*
 * NACA 0012 CT1 case
 */

	psi = -2.89-2.41*sin(*time*2*3.1415926*50.32);
	theta = 0.; 
	phi = 0.;

	printf("Pitch angles is %lf \nRoll angle is %lf \nYaw angle is %lf \n", theta, phi, psi);
	
	tmpfloat[0] = psi;
	tmpfloat[1] = theta;
	tmpfloat[2] = phi;
	
	if(  (TmpNode = m3l_Mklist("RotCenter", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0.25;
	tmpfloat[1] = 0.;
	tmpfloat[2] = 0.;
	
	if(  (TmpNode = m3l_Mklist("TransVec", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0;
	tmpfloat[1] = 0;
	tmpfloat[2] = 0;
	
	
	if(m3l_Cat(Snode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
		Error("CatData");

	client_sender(Snode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);

	
	if(m3l_Umount(&Gnode) != 1)
		Perror("m3l_Umount");
	if(m3l_Umount(&Snode) != 1)
		Perror("m3l_Umount");/* 
 * close socket
 */
	if( close(sockfd) == -1)
		Perror("close");
		

 	}


     return 0; 
}
