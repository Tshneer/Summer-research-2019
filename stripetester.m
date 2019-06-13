clear all
close all

load notes.txt
load lattice_store.csv
load converted_lattice_points.csv


fileID = fopen('notes.txt','r'); formatSpec = '%f';

notes = fscanf(fileID,formatSpec);

length_of_side_of_lattice=notes(1);
collect=notes(2);
scaled_syncronization_parameter=notes(3);
reproductive_rate=notes(4);
radio=notes(5);
dispersal_fraction=notes(6);
test=notes(7);
burn=notes(8);
mean_energy_of_a_patch=notes(9);
number_of_stripes=notes(10);
stripe_test_result=notes(end);


last_snap_shot=lattice2dim(length_of_side_of_lattice,length_of_side_of_lattice,lattice_store);
converted_last_snap_shot=lattice2dim(length_of_side_of_lattice,length_of_side_of_lattice,converted_lattice_points);

ascending_diagonal_stripe_stripe=0;
vertical_stripe=0;
horizontal_stripe=0;
descending_diagonal_stripe=0;
neither_a_striped_nor_a_coherent_state_stripe=0;
coherent_state=0;

if stripe_test_result==1
   ascending_diagonal_stripe_stripe=1;
   disp('ascending_diagonal_strip')
end

if stripe_test_result==2
    horizontal_stripe=1;
     disp('horizontal_strip')
end

if stripe_test_result==3
    vertical_stripe=1;
    disp('vertical_strip')
end

if stripe_test_result==4
    descending_diagonal_stripe=1;
    disp('descending_diagonal_strip')
end

if stripe_test_result==5
    neither_a_striped_nor_a_coherent_state_stripe=1;
    disp('neither_a_striped_nor_a_coherent_state_stripe')
end

if stripe_test_result==6
    coherent_state=1;
    disp('coherent_state')
end

imagesc(last_snap_shot)
%imagesc(converted_last_snap_shot)
daspect([1 1 1])
colorbar('eastoutside')

disp(['scaled_syncronization_parameter is: ',num2str(scaled_syncronization_parameter)])

