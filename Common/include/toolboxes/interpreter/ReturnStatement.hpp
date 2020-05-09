#include "Statement.hpp"

namespace interpreter {

  class ReturnStatement final : public Statement {

  protected:
    calculator expr; //!<\brief Expression of the return statement
    bool value_omitted = true; //!<\brief Boolean indicating that a value has been returned

  protected:

    /*!
     * \brief Compile the statement, i.e. transform the source code into objects.
     * \param code - Pointer to the beginning of the code for the statement
     * \param[out] rest - Pointer to the rest of the code
     * \param parent_scope - Scope of the parent statement
     */
    void _compile(const char* code, const char** rest, TokenMap parent_scope) override;

    /*!
     * \brief Execute the statement
     * \param[in] - The scope of execution
     * \return the return state
     */
    returnState _exec(TokenMap scope) const override;

  public:
    /*!
     * \brief Default constructor
     */
    ReturnStatement() {}

    /*!
     * \brief Constructor using a return type and value
     * \param[in]  code - Pointer to the beginning of the code for the statement
     * \param[in] parent_scope - Scope of the parent statement
     */
    ReturnStatement(const char* code, const char** rest = 0, TokenMap parent_scope = &TokenMap::empty) {
      _compile(code, rest, parent_scope);
    }

    /*!
     * \brief Clone the statement
     * \return pointer to the cloned statement
     */
    Statement* clone() const override {
      return new ReturnStatement(*this);
    }

    /*!
     * \brief literal name of the statement
     */
    static constexpr char literal[] = "return";

  };

}
