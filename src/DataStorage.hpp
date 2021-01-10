#ifndef BIFROST_DATA_STORAGE_HPP
#define BIFROST_DATA_STORAGE_HPP

class DataStorage {

    public:

        DataStorage();
        DataStorage(const vector<string>& color_names_);
        DataStorage(const DataStorage& o);
        DataStorage(DataStorage&& o);

        DataStorage<U>& operator=(const DataStorage& o);
        DataStorage<U>& operator=(DataStorage&& o);

        void clear();

        vector<string> getColorNames(const const_UnitigColorMap<U>& um) const;
        vector<string> getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const;

        bool write(const string& prefix_output_filename, const bool verbose = false) const;
        bool read(const string& filename_colors, const size_t nb_threads = 1, const bool verbose = false);

    private:

        vector<string> color_names;
};

#endif
