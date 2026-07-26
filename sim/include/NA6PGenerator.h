// NA6PCCopyright

#ifndef NA6P_GENERATOR_H
#define NA6P_GENERATOR_H

#include <TObject.h>
#include <fairlogger/Logger.h>
#include <string>
#include <array>

class NA6PMCStack;
class TMethodCall;

class NA6PGenerator : public TObject
{
 public:
  NA6PGenerator(const std::string& name) : mName(name) {}
  NA6PGenerator() = default;
  virtual ~NA6PGenerator() = default;
  virtual void clear() {}
  virtual void setStack(NA6PMCStack* stack) { mStack = stack; }
  virtual void setUserVertexMethod(TMethodCall* m) { mUserVertexMethod = m; }
  auto getUserVertexMethod() const { return mUserVertexMethod; }
  auto getStack() { return mStack; }
  auto getStack() const { return mStack; }

  virtual long canGenerateMaxEvents() const { return -1; } // if no limit, return negative number
  virtual void skipEvents(long) {}
  virtual void generate() = 0;
  virtual void init();
  virtual void generatePrimaryVertex(int maxTrials = 1000000);
  virtual void setOrigin(const std::array<double, 3>& v) { mOrigin = v; }
  template <typename T>
  void setOrigin(T x, T y, T z)
  {
    mOrigin = {x, y, z};
  }
  auto& getOrigin() { return mOrigin; }
  auto& getOrigin() const { return mOrigin; }
  auto getOriginX() const { return mOrigin[0]; }
  auto getOriginY() const { return mOrigin[1]; }
  auto getOriginZ() const { return mOrigin[2]; }

  void setRandomSeed(ULong64_t r) { mRandomSeed = r; }
  auto getRandomSeed() const { return mRandomSeed; }

  bool isInitDone() const { return mInitDone; }
  void setInitDone() { mInitDone = true; }

  const std::string& getName() const { return mName; }

  void setVerbosity(int v) { mVerbosity = v; }
  auto getVerbosity() const { return mVerbosity; }

 protected:
  std::string mName = {};
  std::array<double, 3> mOrigin{}; //! event origin
  NA6PMCStack* mStack = nullptr;   //! externally set stack
  TMethodCall* mUserVertexMethod = nullptr; //! externally set user method for vertex generation
  ULong64_t mRandomSeed = 0;       //  random seed (used if non-0)
  int mVerbosity = 0;              //!
  bool mInitDone = false;

  ClassDef(NA6PGenerator, 3);
};

#endif
