#include "FileSplitter.h"

void FileSplitter::copy(std::istream& is, std::ostream& os, off_t start_offset, size_t end_offset) {
   // naive implementation with one buffer allocation
   size_t size = end_offset - start_offset;
   if (size > buf_size) {
      delete [] buf;
      buf_size = size + 32*1024;
      buf = new char[buf_size];
   }
   is.seekg(start_offset);
   is.read(buf, size);
   os.write(buf, size);
}

void FileSplitter::split(size_t slice_cnt) { // version naive
   size_t instance_info_size = instance_info_v.size();
   if (instance_info_size < slice_cnt) {
      slice_cnt = instance_info_size;
   }
   size_t instance_cnt_per_file = instance_info_size / slice_cnt;
   if ( (instance_info_size % slice_cnt)!=0 ) {
         size_t mintot = instance_cnt_per_file*slice_cnt;
         size_t maxtot = (instance_cnt_per_file+1)*slice_cnt;
	 if ( (instance_info_size - mintot) > (maxtot - instance_info_size) ) {
	    instance_cnt_per_file++;
	 }
   }
   std::ifstream instance_file_is(instance_file.c_str());
   int instpos = instance_file.find(".instance");
   for (size_t nn = 0; nn < slice_cnt; ++nn) {
      std::ostringstream ostr;
      ostr << instance_file.substr(0, instpos) << "_" << nn+1 << ".instance";
      std::string slice_file = ostr.str();
      std::ofstream os(slice_file.c_str());
      // writing header each time
      copy(instance_file_is, os, header_offset, header_size);
      size_t beg = nn*instance_cnt_per_file;
      size_t end = (nn == slice_cnt-1 ? instance_info_size : (nn+1)*instance_cnt_per_file) - 1;

      const InstanceInfo& instance_info_beg = instance_info_v[beg];
      const InstanceInfo& instance_info_end = instance_info_v[end];
      copy(instance_file_is, os, instance_info_beg.offset, instance_info_end.offset + instance_info_end.size);
      os.close();
   }
   // write the effective number of slices in a file
   std::ostringstream ostr2;
   ostr2 << instance_file.substr(0, instpos) << ".slicecount";
   std::string slice2_file = ostr2.str();
   std::ofstream os2(slice2_file.c_str()); 
   os2<<"@Total Number of Effective Slices\n"<<slice_cnt<<"\t";
   os2.close();
}

// LATER
/*
   void FileSplitter::split2(size_t slice_cnt) {
   InstanceInfo& g_instance_info_beg = instance_info_v[0];
   InstanceInfo& g_instance_info_end = instance_info_v[instance_info_v.size()-1];
   size_t all_instance_size = g_instance_info_end.offset + g_instance_info_end.size - g_instance_info_beg.offset;
   size_t instance_info_size = instance_info_v.size();
   size_t size_per_slice = all_instance_size / slice_cnt;
   int ind_beg = 0;
   size_t size_sum = 0;
   int slice_num = 0;
   for (int ind = 0; ind < instance_info_size; ++ind) {
   const InstanceInfo& instance_info_cur = instance_info_v[ind];
   if (ind == instance_info_size - 1 || size_sum + instance_info_cur.size >= size_per_splice) {
   std::string slice_file = instance_file + "_" + slice_num;
   std::oftream os(slice_file);
   const InstanceInfo& instance_info_beg = instance_info_v[ind_beg];
   copy(instance_file_is, os, instance_info_beg.offset, instance_info_cur.offset + instance_info_cur.size);
   ind_beg = ind+1;
   slice_num++;
   size_sum = 0;
   } else {
   size_sum += instance_info_cur.size;
   }
   }
   }
   */
