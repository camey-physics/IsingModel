
class Model {
public:
    virtual ~Model() = default;
    virtual void initializeState() = 0;
    virtual void copyStateFrom(const Model& other) = 0;  // Pure virtual

    virtual double calcEnergy() const = 0;
    virtual void updateSweep(int numSweeps) = 0;
    // not sure about this yet; there may be a better way to create clones if I keep it 
    // virtual std::unique_ptr<Model> clone() const = 0; 

};