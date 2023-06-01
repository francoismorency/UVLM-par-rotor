%% Save To Excel

filename0='Thrust_Data.csv';
filename1='ReY.csv';
filename2='CLY_Values.csv';
filename3='CLYJ_Values.csv';
filename4='Initial_Data.csv';
filename5='Wake_Shape_x.csv';
filename6='Wake_Shape_y.csv';
filename7='Wake_Shape_z.csv';
filename8='Wake_ShapeA_x.csv';
filename9='Wake_ShapeA_y.csv';
filename10='Wake_ShapeA_z.csv';

DD=[Azimuth];DD1=[NREV];DD2=[CT];DD3=[CQ];DD4=[CL];DD5=[FM];DD6=[IterationTime];
Data=[DD;DD1;DD2;DD3;DD4;DD5;DD6];
Initial_Data=[C B RootY IB JB RPM ALFA*180/pi AirfoilShape DTETA NREVSLOWSTART SpeedofSound Solidity CH ClimbRatio];

csvwrite(filename0,Data);
csvwrite(filename1,ReY);
csvwrite(filename2,CLY);
csvwrite(filename3,CLYJ);
csvwrite(filename4,Initial_Data);
csvwrite(filename5,QW(:,:,1));
csvwrite(filename6,QW(:,:,2));
csvwrite(filename7,QW(:,:,3));

csvwrite(filename8,QWA(:,:,1));
csvwrite(filename9,QWA(:,:,2));
csvwrite(filename10,QWA(:,:,3));
