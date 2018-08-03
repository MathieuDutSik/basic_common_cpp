#ifndef GMP_SERIALIZATION_INCLUDE
#define GMP_SERIALIZATION_INCLUDE


namespace boost { namespace serialization {

    // mpq_class
    
    template<class Archive>
      inline void load(Archive & ar,
		       mpq_class & val,
		       const unsigned int version)
    {
      std::cerr << "load(mpq_class), step 1\n";
      std::string str;
      ar & make_nvp("mpq", str);
      std::istringstream is(str);
      is >> val;
      std::cerr << "load(mpq_class), step 2\n";
    }


    
    template<class Archive>
      inline void save(Archive & ar,
		       mpq_class const& val,
		       const unsigned int version)
    {
      std::cerr << "save(mpq_class), step 1\n";
      std::ostringstream os;
      os << val;
      std::string str=os.str();
      ar & make_nvp("mpq", str);
      std::cerr << "save(mpq_class), step 2\n";
    }

    template<class Archive>
      inline void serialize(Archive & ar,
			    mpq_class & val,
			    const unsigned int version)
    {
      std::cerr << "split_free(mpq_class), step 1\n";
      split_free(ar, val, version);
      std::cerr << "split_free(mpq_class), step 2\n";
    }

    // mpz_class

    template<class Archive>
      inline void load(Archive & ar,
		       mpz_class & val,
		       const unsigned int version)
    {
      std::cerr << "load(mpz_class), step 1\n";
      std::string str;
      ar & make_nvp("mpz", str);
      std::istringstream is(str);
      is >> val;
      std::cerr << "load(mpz_class), step 2\n";
    }


    
    template<class Archive>
      inline void save(Archive & ar,
		       mpz_class const& val,
		       const unsigned int version)
    {
      std::cerr << "save(mpz_class), step 1\n";
      std::ostringstream os;
      os << val;
      std::string str=os.str();
      ar & make_nvp("mpz", str);
      std::cerr << "save(mpz_class), step 2\n";
    }

    template<class Archive>
      inline void serialize(Archive & ar,
			    mpz_class & val,
			    const unsigned int version)
    {
      std::cerr << "split_free(mpz_class), step 1\n";
      split_free(ar, val, version);
      std::cerr << "split_free(mpz_class), step 2\n";
    }
}}



#endif
