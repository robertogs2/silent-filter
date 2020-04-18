/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <exception>

#ifndef SILENT_EXCEPTION_HPP
#define SILENT_EXCEPTION_HPP

namespace silent {

  /**
   * Exception
   */
  class Exception : public std::exception {
    std::string _name;
  public:
    inline Exception(const std::string& name="silent exception") : _name(name) { }
    inline const char* what() const noexcept { return _name.c_str(); }
  };
  
}

#endif

