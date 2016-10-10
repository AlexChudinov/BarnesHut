#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <utility>
#include <stack>

namespace data_structs {

using std::pair;
using std::stack;
template<class data, class branch>
using triple = pair<data, pair<branch, branch>>;

template<class T> class binary_tree;

/**
 * Simple binary tree impementation
 */
template<class T> class binary_tree_node
{
    typedef binary_tree_node<T> node;

    T data;
    node* leftSon;
    node* rightSon;
    friend class binary_tree<T>;

    /**
     * default constructor
     */
    binary_tree_node() : leftSon(nullptr), rightSon(nullptr) {}

    /**
     * Creates a new node decomposing a type T1 value
     */
    template<class FunObj, class T1>
    binary_tree_node(const FunObj& fun, const T1& t)
        :
          leftSon(nullptr), rightSon(nullptr)
    {
        typedef pair<T1, node*> pair;

        //Push first vector to a stack
        stack<pair> dfsStack;
        dfsStack.push(pair(t,this));

        while(!dfsStack.empty())
        {
            T1 t = dfsStack.top().first;
            node * pnode = dfsStack.top().second;
            dfsStack.pop();

            //Split inputed data and returns a node value
            triple<T, T1> nodeData = fun(t);
            pnode->data = nodeData.first;

            if(!nodeData.second.first.empty())
            {
                this->leftSon = new node();
                dfsStack.push(pair(nodeData.second.first,
                                   this->leftSon));
            }
            if(!nodeData.second.second.empty())
            {
                this->rightSon = new node();
                dfsStack.push(pair(nodeData.second.second,
                                   this->rightSon));
            }
        }
    }

    ~binary_tree_node()
    {
        stack<node> dfsStack;
        dfsStack.push(this);

        while(!dfsStack.empty())
        {
            node* pNode = dfsStack.top();
            dfsStack.pop();

            if (pNode->leftSon) dfsStack.push(pNode->leftSon);
            if (pNode->rightSon)dfsStack.push(pNode->rightSon);

            delete pNode;
        }
    }
};

template<class T> class binary_tree
{
    typedef binary_tree_node<T> node;

    node* root;
public:

    template<class FunObj, class T1>
    binary_tree(const FunObj& fun, const T1& t)
        :
          root(new node(fun, t)){}

    void clear() { delete root; }
};

}

#endif // BINARYTREE_H
