LUTs not uploaded as there were too many.
	
Changing distribution:
In Report things 2 Random length and doping
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts
	○ exp_length_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts
	○ poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts
	○ "C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts.raw"
		
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2
	○ exp_length_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2
	○ poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3 (LUT: Test3)
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 4 (LUT: Test4)
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 5 (LUT: Test5)
	○ poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 5.txt"
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 6 (LUT: Test5)
 ... up to 11
 (Stopped using different LUT prefixes as it allowed them to be reused)

Changing inputs and outputs:
In Report things 2 Random length and doping
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts MIMG (0-4,5-9,10-15) (14) (LUT: Test2)
	○ exp_length_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts
	○ poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts
	○ "left":([range(4), range(5,9), range(10,15)],["V1","V2","V3"])
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG (0-4,5-9,10-15) (14) 3 (LUT: Test3)
	○ exp_length_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3
	○ poisson_doping_dist_Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3


- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts MIMG (0-5,10-15) (0,14)

- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 1 (test 2) 
(This was to see what happened if you changed 5 inputs along one edge eg. 10000 to 00011, and one output at node 14, see below for input location)
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 2
	○ "left":([range(3), range(3,6), range(6,9), range(9,12), range(12,15)],["V1","V2","V3","V4","V5"])
- Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 5 (LUT:Test5)
	
