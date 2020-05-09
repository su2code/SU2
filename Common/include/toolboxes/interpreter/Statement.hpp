#include "shunting-yard.h"
#include "shunting-yard-exceptions.h"
#pragma once

namespace interpreter {

  extern GlobalScope globalScope; //!<\brief The global function scope

  /*!
   * \brief The return types a statement can have
   */
  enum returnType {
    NORMAL,
    FINISH,
    RETURN,
    YIELD,
    BREAK,
    CONTINUE,
    THROW
  };

  /*!
   * \brief A structure to represent the type and value that is returned by a statement
   */
  struct returnState {
    returnType type;           //!<\brief The return type
    packToken value;           //!<\brief The value that is returned
    bool value_omitted = true; //!<\brief Boolean indicating whether a value has been returned

    /*!
     * \brief Default constructor
     */
    returnState() : type(NORMAL), value(packToken::None()) {}

    /*!
     * \brief Constructor using a return type
     * \param[in] - The return type
     */
    returnState(const returnType& type) : type(type), value(packToken::None()) {}

    /*!
     * \brief Constructor using a return value
     * \param[in] value - The value
     */
    returnState(packToken value) : type(NORMAL), value(value) {}

    /*!
     * \brief Constructor using a return type and value
     * \param[in] type - The return type
     * \param[in] value - The return value
     * \param[in] value_omitted - Boolean indicating whether a value has been returned
     */
    returnState(const returnType& type, packToken value, bool value_omitted = true)
      : type(type), value(value), value_omitted(value_omitted) {}
  };

  /*!
   * \brief Abstract class to represent a statement
   */
  class Statement {
  protected:
    /*!
     * \brief Compile the statement, i.e. transform the source code into objects.
     * \param[in] code - Pointer to the beginning of the code for the statement
     * \param[out] rest - Pointer to the rest of the code
     * \param[in] parent_scope - Scope of the parent statement
     */
    virtual void _compile(const char* code, const char** rest, TokenMap parent_scope) = 0;

    /*!
     * \brief Execute the statement
     * \param[in] - The scope of execution
     * \return the return state
     */
    virtual returnState _exec(TokenMap scope) const = 0;

  public:

    /*!
     * \brief Virtual desctructor
     */
    virtual ~Statement() {}

    /*!
     * \brief Compile the statement, i.e. transform the source code into objects.
     * \param[in]  code - Pointer to the beginning of the code for the statement
     * \param[out] rest - Pointer to the rest of the code
     * \param[in]  parent_scope - Scope of the parent statement
     */
    void compile(const char* code, const char** rest = 0, TokenMap parent_scope = &TokenMap::empty) {
      return _compile(code, rest, parent_scope);
    }

    /*!
     * \brief Execute the statement
     * \param[in] - The scope of execution
     * \return the return state
     */
    returnState exec(TokenMap scope) const { return _exec(scope); }

    /*!
     * \brief Clone the statement
     * \return pointer to the cloned statement
     */
    virtual Statement* clone() const = 0;

    /*!
     * \brief Build a statement from source code
     * \param source - Pointer to the source code
     * \param scope - Scope where the statement is embedded
     * \return - Pointer to the created statement
     */
    static Statement* buildStatement(const char** source, TokenMap scope);

  };

}
