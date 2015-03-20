#include <sys/types.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>


class FileSplitter {
   off_t header_offset;
   size_t header_size;

   struct InstanceInfo {
      off_t offset;
      size_t size;
      // stats eventuelles: # de types etc.
      InstanceInfo() : offset(0), size(0) {}
      InstanceInfo(off_t offset) : offset(offset), size(0) {}
   };

   std::vector<InstanceInfo> instance_info_v; // vecteur d'objets de type InstanceInfo
   std::string instance_file;
   std::ostream& instance_os;
   size_t buf_size;
   char* buf;

   public:

   FileSplitter(const std::string& instance_file, std::ostream& instance_os) : instance_file(instance_file), instance_os(instance_os) {
      header_offset = (off_t)-1;
      header_size = (off_t)-1;
      buf_size = 0;
      buf = NULL;
   }

   void startWritingHeader() {
      header_offset = instance_os.tellp();
   }

   void endWritingHeader() {
      header_size = instance_os.tellp() - header_offset;
   }

   void startWritingInstance() {
      InstanceInfo instance_info(instance_os.tellp());
      instance_info_v.push_back(instance_info);
   }

   void endWritingInstance(/* ... */) {
      InstanceInfo& instance_info = instance_info_v[instance_info_v.size()-1];
      instance_info.size = instance_os.tellp() - instance_info.offset; // attention piquet !!!
      // instance_info.stats = ...
   }

   void copy(std::istream& is, std::ostream& os, off_t start_offset, size_t end_offset);

   void split(size_t slice_cnt);

};

