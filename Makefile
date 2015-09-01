SRC_DIR = ./src
SRC__f_DIR = ./src_FastPCAdapt
SRC__c_DIR = ./src_convert
OBJ_DIR = ./obj
OBJ__c_DIR = ./obj_convert
lapack_SRC_DIR = ./src_Lapack
lapack_OBJ_DIR = ./obj_Lapack
PROG = PCAdapt  
PROG_c1 = vcf2pcadapt
PROG_c2 = ped2pcadapt
OBJS = $(OBJ_DIR)/Data.o $(OBJ_DIR)/random.o $(OBJ_DIR)/mcmc.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/updates.o $(OBJ_DIR)/linAlgebra.o $(OBJ_DIR)/stat.o $(OBJ_DIR)/FastPCAdapt.o $(OBJ_DIR)/Cov_line.o $(OBJ_DIR)/Load_line.o $(OBJ_DIR)/Data__f.o  $(lapack_OBJ_DIR)/*.o 
OBJS_c = $(OBJ__c_DIR)/ancestrymap.o $(OBJ__c_DIR)/error_matrix.o $(OBJ__c_DIR)/geno.o $(OBJ__c_DIR)/io_data_double.o $(OBJ__c_DIR)/io_data_float.o $(OBJ__c_DIR)/io_data_int.o $(OBJ__c_DIR)/io_error.o $(OBJ__c_DIR)/io_tools.o $(OBJ__c_DIR)/ped.o $(OBJ__c_DIR)/print_bar.o $(OBJ__c_DIR)/register_convert.o $(OBJ__c_DIR)/vcf2geno.o
CC = gcc
CFLAGS = -g -lm -O3 -Wformat=0 -Wno-div-by-zero 

$(PROG): $(OBJS) $(OBJS_c)
	$(CC) $(SRC_DIR)/PCAdapt.c -o $@ $(OBJS) $(CFLAGS)
	$(CC) $(SRC__c_DIR)/main_vcf2geno.c -o $(PROG_c1) $(OBJS_c) $(CFLAGS)
	$(CC) $(SRC__c_DIR)/main_ped2geno.c -o $(PROG_c2) $(OBJS_c) $(CFLAGS)

$(OBJ_DIR)/Data.o: $(SRC_DIR)/Data.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/Data.c -o $(OBJ_DIR)/Data.o

$(OBJ_DIR)/mcmc.o: $(SRC_DIR)/mcmc.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/mcmc.c -o $(OBJ_DIR)/mcmc.o

$(OBJ_DIR)/random.o: $(SRC_DIR)/random.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/random.c -o $(OBJ_DIR)/random.o

$(OBJ_DIR)/matrix.o: $(SRC_DIR)/matrix.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/matrix.c -o $(OBJ_DIR)/matrix.o

$(OBJ_DIR)/linAlgebra.o: $(SRC_DIR)/linAlgebra.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/linAlgebra.c -o $(OBJ_DIR)/linAlgebra.o

$(OBJ_DIR)/updates.o: $(SRC_DIR)/updates.h
	$(CC) $(CFLAGS) -c $(SRC_DIR)/updates.c -o $(OBJ_DIR)/updates.o

$(OBJ_DIR)/FastPCAdapt.o: $(SRC__f_DIR)/FastPCAdapt.h
	$(CC) $(CFLAGS) -c $(SRC__f_DIR)/FastPCAdapt.c -o $(OBJ_DIR)/FastPCAdapt.o

$(OBJ_DIR)/Data__f.o: $(SRC__f_DIR)/Data__f.h
	$(CC) $(CFLAGS) -c $(SRC__f_DIR)/Data__f.c -o $(OBJ_DIR)/Data__f.o

$(OBJ_DIR)/Cov_line.o: $(SRC__f_DIR)/Cov_line.h
	$(CC) $(CFLAGS) -c $(SRC__f_DIR)/Cov_line.c -o $(OBJ_DIR)/Cov_line.o

$(OBJ_DIR)/Load_line.o: $(SRC__f_DIR)/Load_line.h
	$(CC) $(CFLAGS) -c $(SRC__f_DIR)/Load_line.c -o $(OBJ_DIR)/Load_line.o

$(OBJ_DIR)/stat.o: $(SRC__f_DIR)/stat.h
	$(CC) $(CFLAGS) -c $(SRC__f_DIR)/stat.c -o $(OBJ_DIR)/stat.o

$(OBJ__c_DIR)/ancestrymap.o: $(SRC__c_DIR)/ancestrymap.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/ancestrymap.c -o $(OBJ__c_DIR)/ancestrymap.o

$(OBJ__c_DIR)/error_matrix.o: $(SRC__c_DIR)/error_matrix.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/error_matrix.c -o $(OBJ__c_DIR)/error_matrix.o

$(OBJ__c_DIR)/geno.o: $(SRC__c_DIR)/geno.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/geno.c -o $(OBJ__c_DIR)/geno.o

$(OBJ__c_DIR)/io_data_double.o: $(SRC__c_DIR)/io_data_double.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/io_data_double.c -o $(OBJ__c_DIR)/io_data_double.o

$(OBJ__c_DIR)/io_data_float.o: $(SRC__c_DIR)/io_data_float.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/io_data_float.c -o $(OBJ__c_DIR)/io_data_float.o

$(OBJ__c_DIR)/io_data_int.o: $(SRC__c_DIR)/io_data_int.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/io_data_int.c -o $(OBJ__c_DIR)/io_data_int.o

$(OBJ__c_DIR)/io_error.o: $(SRC__c_DIR)/io_error.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/io_error.c -o $(OBJ__c_DIR)/io_error.o

$(OBJ__c_DIR)/io_tools.o: $(SRC__c_DIR)/io_tools.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/io_tools.c -o $(OBJ__c_DIR)/io_tools.o

$(OBJ__c_DIR)/ped.o: $(SRC__c_DIR)/ped.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/ped.c -o $(OBJ__c_DIR)/ped.o

$(OBJ__c_DIR)/print_bar.o: $(SRC__c_DIR)/print_bar.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/print_bar.c -o $(OBJ__c_DIR)/print_bar.o

$(OBJ__c_DIR)/register_convert.o: $(SRC__c_DIR)/register_convert.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/register_convert.c -o $(OBJ__c_DIR)/register_convert.o

$(OBJ__c_DIR)/vcf2geno.o: $(SRC__c_DIR)/vcf2geno.h
	$(CC) $(CFLAGS) -c $(SRC__c_DIR)/vcf2geno.c -o $(OBJ__c_DIR)/vcf2geno.o

lapack: 
	$(CC) $(CFLAGS) -c $(lapack_SRC_DIR)/*.c
	mv *.o $(lapack_OBJ_DIR)/

clean:
	rm -f $(OBJ_DIR)/*.o $(PROG) 
	rm -f $(OBJ__c_DIR)/*.o $(PROG_c1) $(PROG_c2)

realclean:
	rm -f $(OBJ_DIR)/*.o $(PROG) 
	rm -f $(OBJ__c_DIR)/*.o $(PROG_c1) $(PROG_c2)
	rm -f $(lapack_OBJ_DIR)/*.o $(OBJ_DIR)/*.o $(PROG) 

