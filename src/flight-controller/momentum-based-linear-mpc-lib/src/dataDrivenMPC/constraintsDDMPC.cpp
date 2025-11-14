#include <FlightControlUtils.h>
#include <dataDrivenMPC/DDconstant.h>
#include <dataDrivenMPC/constraintsDDMPC.h>

namespace DDMPC
{

ConstraintSystemDynamicDD::ConstraintSystemDynamicDD(
    const int nVar, const int nStates, const int nJoints, const int nThrottle, const int nIter)
    : IQPConstraint(nVar, nStates * nIter)
{
    m_systemDynamicDF = std::make_unique<SystemDynamicDD>(nStates, nJoints, nThrottle);
    m_A.resize(nStates, nStates);
    m_BJoints.resize(nStates, nJoints);
    m_BThrottle.resize(nStates, nThrottle);
    m_c.resize(nStates);
    m_nStates = nStates;
    m_nJoints = nJoints;
    m_nThrottle = nThrottle;
    m_nIter = nIter;
}

const bool ConstraintSystemDynamicDD::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("nIterSmall", m_nSmallSteps))
    {
        yError() << "Parameter 'nIterSmall' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("controlHorizon", m_ctrlHorizon))
    {
        yError() << "Parameter 'controlHorizon' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("periodMPCSmallSteps", m_deltaTSmallSteps))
    {
        yError() << "Parameter 'periodMPCSmallSteps' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("periodMPCLargeSteps", m_deltaTLargeSteps))
    {
        yError() << "Parameter 'periodMPCLargeSteps' not found in the config file.";
        return false;
    }
    // the warp function is computed as: w(tau) = beta1 * tau + beta2 * tau^2
    // the coefficients beta1 and beta2 are computed solving the linear system:
    // w(1) - w(0) = deltaTSmallSteps
    // w(nSmallSteps) = 0.1
    m_beta2 = (m_deltaTLargeSteps - m_nSmallSteps * m_deltaTSmallSteps)
              / (m_nSmallSteps * (m_nSmallSteps - 1));
    m_beta1 = m_deltaTSmallSteps - m_beta2;
    return m_systemDynamicDF->configure(parametersHandler, qpInput);
}

void ConstraintSystemDynamicDD::configureDynVectorsSize(QPInput& qpInput)
{
    m_updateThrottle = true;
    m_counter = 0;
}

const bool ConstraintSystemDynamicDD::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    if (!m_systemDynamicDF->updateDynamicMatrices(qpInput))
    {
        yError() << "ConstraintSystemDynamicDD::computeConstraintsMatrixAndBounds: error in "
                    "updateDynamicMatrices";
        return false;
    }
    m_systemDynamicDF->getAMatrix(m_A);
    m_systemDynamicDF->getBJointsMatrix(m_BJoints);
    m_systemDynamicDF->getBThrottleMatrix(m_BThrottle);
    m_systemDynamicDF->getCVector(m_c);
    m_linearMatrix.setZero();
    m_lowerBound.setZero();
    m_upperBound.setZero();
    for (int i = 0; i < m_nIter; i++)
    {
        if (i < m_nSmallSteps)
        {
            m_deltaT = warp_function(i + 1) - warp_function(i);
        } else
        {
            m_deltaT = m_deltaTLargeSteps;
        }
        m_linearMatrix.block(i * m_nStates, i * m_nStates, m_nStates, m_nStates)
            = Eigen::MatrixXd::Identity(m_nStates, m_nStates) + m_deltaT * m_A;
        m_linearMatrix.block(i * m_nStates, (i + 1) * m_nStates, m_nStates, m_nStates)
            = -Eigen::MatrixXd::Identity(m_nStates, m_nStates);
        if (i < m_ctrlHorizon)
        {
            m_linearMatrix.block(i * m_nStates,
                                 m_nStates * (m_nIter + 1) + i * m_nJoints,
                                 m_nStates,
                                 m_nJoints)
                = m_deltaT * m_BJoints;
        } else
        {
            m_linearMatrix.block(i * m_nStates,
                                 m_nStates * (m_nIter + 1) + (m_ctrlHorizon - 1) * m_nJoints,
                                 m_nStates,
                                 m_nJoints)
                = m_deltaT * m_BJoints;
        }
        if (i < m_nSmallSteps)
        {
            m_linearMatrix.block(i * m_nStates,
                                 m_nStates * (m_nIter + 1) + m_ctrlHorizon * m_nJoints,
                                 m_nStates,
                                 m_nThrottle)
                = m_deltaT * m_BThrottle;

        } else if (i >= m_nSmallSteps && i < m_ctrlHorizon)
        {
            m_linearMatrix.block(i * m_nStates,
                                 m_nStates * (m_nIter + 1) + m_ctrlHorizon * m_nJoints
                                     + (i - (m_nSmallSteps - 1)) * m_nThrottle,
                                 m_nStates,
                                 m_nThrottle)
                = m_deltaT * m_BThrottle;
        } else
        {
            m_linearMatrix.block(i * m_nStates,
                                 m_nStates * (m_nIter + 1) + m_ctrlHorizon * m_nJoints
                                     + (m_ctrlHorizon - (m_nSmallSteps)) * m_nThrottle,
                                 m_nStates,
                                 m_nThrottle)
                = m_deltaT * m_BThrottle;
        }
        m_lowerBound.segment(i * m_nStates, m_nStates) = -m_deltaT * m_c;
        m_upperBound.segment(i * m_nStates, m_nStates) = -m_deltaT * m_c;
    }
    if (m_counter == m_nSmallSteps)
    {
        m_updateThrottle = true;
        m_counter = 0;
    } else
    {
        m_updateThrottle = false;
    }
    m_counter++;
    return true;
}

const bool
ConstraintSystemDynamicDD::populateVectorsCollection(QPInput& qpInput,
                                                     const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

const bool ConstraintSystemDynamicDD::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

double ConstraintSystemDynamicDD::warp_function(double t)
{
    return m_beta1 * t + m_beta2 * t * t;
}

ConstraintInitialStateDD::ConstraintInitialStateDD(const unsigned int nStates,
                                                   const unsigned int nVar)
    : IQPConstraintInitialState(nStates, nVar)
{
}

const bool ConstraintInitialStateDD::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("periodMPC", m_periodMPC))
    {
        yError() << "Parameter 'periodMPC' not found in the config file.";
        return false;
    }
    return true;
}

void ConstraintInitialStateDD::configureDynVectorsSize(QPInput& qpInput)
{
    m_robot = qpInput.getRobot();
    m_vectorsCollectionServer = qpInput.getVectorsCollectionServer();
    m_initialThrust.resize(m_robot->getNJets());
    m_initialThrustDot.resize(m_robot->getNJets());
    m_initialThrust.resize(m_robot->getNJets());
    m_initialThrustDot.resize(m_robot->getNJets());
    m_initialThrust = m_robot->getJetThrusts();
    m_initialThrustDot.setZero();
    m_initialAngMom = m_robot->getMomentum(true).bottomRows(3);
    m_initialCoMPos = m_robot->getPositionCoM();
    m_initialLinMom = m_robot->getMomentum(true).topRows(3);
    m_initialRPY = iDynTree::toEigen(m_robot->getBasePose().getRotation().asRPY());
    m_rpyOld = m_initialRPY;
    m_RPYinit = m_initialRPY;
    m_nTurns.setZero();
    m_initialAlphaGravity = 0.08;
}

bool ConstraintInitialStateDD::updateInitialState(QPInput& qpInput)
{
    this->unwrapRPY();
    m_initialState.setZero();
    m_initialState.segment(CoMPosIdx[0], CoMPosIdx.size()) = m_robot->getPositionCoM();
    m_initialState.segment(linMomIdx[0], linMomIdx.size()) = m_robot->getMomentum(true).topRows(3);
    m_initialState.segment(rpyIdx[0], rpyIdx.size()) = m_rpyUnwrapped;
    m_initialState.segment(angMomIdx[0], angMomIdx.size())
        = m_robot->getMomentum(true).bottomRows(3);
    return true;
}

void ConstraintInitialStateDD::unwrapRPY()
{
    for (int i = 0; i < 3; i++)
    {
        if (m_robot->getBasePose().getRotation().asRPY()(i) - m_rpyOld(i) > M_PI)
        {
            m_nTurns(i)--;
        } else if (m_robot->getBasePose().getRotation().asRPY()(i) - m_rpyOld(i) < -M_PI)
        {
            m_nTurns(i)++;
        }
    }
    m_rpyUnwrapped
        = iDynTree::toEigen(m_robot->getBasePose().getRotation().asRPY()) + 2 * M_PI * m_nTurns;
    m_rpyOld = iDynTree::toEigen(m_robot->getBasePose().getRotation().asRPY());
}

const bool ConstraintInitialStateDD::configureVectorsCollectionServer(QPInput& qpInput)
{
    m_vectorsCollectionServer->populateMetadata("measures::unwrappedRPY", {"roll", "pitch", "yaw"});
    return true;
}

const bool
ConstraintInitialStateDD::populateVectorsCollection(QPInput& qpInput,
                                                    const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    m_vectorsCollectionServer->populateData("measures::unwrappedRPY", m_rpyUnwrapped);
    return true;
}

ThrottleConstraintDD::ThrottleConstraintDD(const int nVar,
                                           const int nStates,
                                           const int nIter,
                                           const int nSmallSteps,
                                           const int ctrlHorizon,
                                           const int throttleStartIdx)
    : IQPConstraint(nVar, N_THRUSTS * (ctrlHorizon - nSmallSteps + 1))
{
    m_nIter = nIter;
    m_nSmallSteps = nSmallSteps;
    m_nStates = nStates;
    m_ctrlHorizon = ctrlHorizon;
    m_nThrottleIdxStart = throttleStartIdx;
}

const bool ThrottleConstraintDD::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("throttleMax", m_throttleMaxValue))
    {
        yError() << "Parameter 'throttleMax' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("throttleMin", m_throttleMinValue))
    {
        yError() << "Parameter 'throttleMin' not found in the config file.";
        return false;
    }
    double periodMPCLargeSteps;
    if (!ptr->getParameter("periodMPCLargeSteps", periodMPCLargeSteps))
    {
        yError() << "Parameter 'periodMPCLargeSteps' not found in the config file.";
        return false;
    }
    double periodMPCSmallSteps;
    if (!ptr->getParameter("periodMPCSmallSteps", periodMPCSmallSteps))
    {
        yError() << "Parameter 'periodMPCSmallSteps' not found in the config file.";
        return false;
    }
    m_ratioSmallLargeStepsPeriod = round(periodMPCLargeSteps / periodMPCSmallSteps);
    return true;
}

void ThrottleConstraintDD::configureDynVectorsSize(QPInput& qpInput)
{
    m_jetModel = qpInput.getJetModel();
    m_throttleMaxValue = m_jetModel->compute_u2v(m_throttleMaxValue);
    m_throttleMinValue = m_jetModel->compute_u2v(m_throttleMinValue);
    m_nJoints = deltaJointIdx.size();
    m_nThrottle = throttleIdx.size();
    m_counter = m_ratioSmallLargeStepsPeriod - 1;
}

const bool ThrottleConstraintDD::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    m_linearMatrix.setZero();
    m_lowerBound.setZero();
    m_upperBound.setZero();
    for (int i = 0; i < m_ctrlHorizon - m_nSmallSteps + 1; i++)
    {
        m_linearMatrix
            .block(i * m_nThrottle, m_nThrottleIdxStart + i * m_nThrottle, m_nThrottle, m_nThrottle)
            .setIdentity();
        if (m_counter != (m_ratioSmallLargeStepsPeriod - 1) && i == 0)
        {
            for (int j = 0; j < N_THRUSTS; j++)
            {
                m_lowerBound[j] = m_jetModel->compute_u2v(qpInput.getThrottleMPC()[j]);
                m_upperBound[j] = m_jetModel->compute_u2v(qpInput.getThrottleMPC()[j]);
            }
        } else
        {
            m_lowerBound.segment(i * m_nThrottle, m_nThrottle).setConstant(m_throttleMinValue);
            m_upperBound.segment(i * m_nThrottle, m_nThrottle).setConstant(m_throttleMaxValue);
        }
    }
    if (m_counter == (m_ratioSmallLargeStepsPeriod - 1))
    {
        qpInput.setUpdateThrottle(true);
        m_counter = 0;
    } else
    {
        qpInput.setUpdateThrottle(false);
        m_counter++;
    }
    return true;
}

const bool ThrottleConstraintDD::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
ThrottleConstraintDD::populateVectorsCollection(QPInput& qpInput,
                                                const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

JointPositionConstraintDD::JointPositionConstraintDD(
    const int nVar, const int nStates, const int nJoints, const int nIter, const int controlHorizon)
    : IQPConstraint(nVar, nJoints * controlHorizon)
{
    m_nIter = nIter;
    m_nStates = nStates;
    m_nJoints = nJoints;
}

const bool JointPositionConstraintDD::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    m_jointPositionMax.resize(m_nJoints);
    m_jointPositionMin.resize(m_nJoints);
    if (!ptr->getParameter("jointPos_max", m_jointPositionMax))
    {
        yError() << "Parameter 'jointPos_max' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("jointPos_min", m_jointPositionMin))
    {
        yError() << "Parameter 'jointPos_min' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("controlHorizon", m_ctrlHorizon))
    {
        yError() << "Parameter 'controlHorizon' not found in the config file.";
        return false;
    }
    // convert the joint position limits from degrees to radians
    m_jointPositionMax = m_jointPositionMax * M_PI / 180.0;
    m_jointPositionMin = m_jointPositionMin * M_PI / 180.0;
    return true;
}

void JointPositionConstraintDD::configureDynVectorsSize(QPInput& qpInput)
{
    m_robot = qpInput.getRobot();
    m_linearMatrix.setZero();
    m_firstIteriation = true;
}

const bool JointPositionConstraintDD::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    m_lowerBound.setZero();
    m_upperBound.setZero();
    for (int i = 0; i < m_ctrlHorizon; i++)
    {
        m_linearMatrix
            .block(i * m_nJoints, m_nStates * (m_nIter + 1) + i * m_nJoints, m_nJoints, m_nJoints)
            .setIdentity();
        m_lowerBound.segment(i * m_nJoints, m_nJoints)
            = m_jointPositionMin - qpInput.getOutputQPJointsPosition().segment(3, m_nJoints);
        m_upperBound.segment(i * m_nJoints, m_nJoints)
            = m_jointPositionMax - qpInput.getOutputQPJointsPosition().segment(3, m_nJoints);
    }
    return true;
}

const bool JointPositionConstraintDD::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
JointPositionConstraintDD::populateVectorsCollection(QPInput& qpInput,
                                                     const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

ThrustContraintDD::ThrustContraintDD(const int nVar,
                                     const int nStates,
                                     const int nIter,
                                     const int nSmallSteps,
                                     const int ctrlHorizon,
                                     const int thrustInitPosition)
    : IQPConstraint(nVar, N_THRUSTS * (ctrlHorizon - nSmallSteps + 1))
{
    m_nIter = nIter;
    m_nSmallSteps = nSmallSteps;
    m_nStates = nStates;
    m_ctrlHorizon = ctrlHorizon;
    m_thrustInitPosition = thrustInitPosition;
}

const bool ThrustContraintDD::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("thrustMax", m_thrustMaxValue))
    {
        yError() << "Parameter 'thrustMax' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("thrustMin", m_thrustMinValue))
    {
        yError() << "Parameter 'thrustMin' not found in the config file.";
        return false;
    }
    double periodMPCLargeSteps;
    if (!ptr->getParameter("periodMPCLargeSteps", periodMPCLargeSteps))
    {
        yError() << "Parameter 'periodMPCLargeSteps' not found in the config file.";
        return false;
    }
    double periodMPCSmallSteps;
    if (!ptr->getParameter("periodMPCSmallSteps", periodMPCSmallSteps))
    {
        yError() << "Parameter 'periodMPCSmallSteps' not found in the config file.";
        return false;
    }
    m_ratioSmallLargeStepsPeriod = round(periodMPCLargeSteps / periodMPCSmallSteps);
    std::cout << "ratioSmallLargeStepsPeriod: " << m_ratioSmallLargeStepsPeriod << std::endl;
    return true;
}
void ThrustContraintDD::configureDynVectorsSize(QPInput& qpInput)
{
    m_thrustMaxValue = (m_thrustMaxValue);
    m_thrustMinValue = (m_thrustMinValue);
    m_robot = qpInput.getRobot();
    m_nJoints = deltaJointIdx.size();
    m_nThrottle = throttleIdx.size();
    m_counter = m_ratioSmallLargeStepsPeriod - 1;
}

const bool ThrustContraintDD::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    m_linearMatrix.setZero();
    m_lowerBound.setConstant(m_thrustMinValue);
    m_upperBound.setConstant(m_thrustMaxValue);
    for (int i = 0; i < m_ctrlHorizon - m_nSmallSteps + 1; i++)
    {
        m_linearMatrix
            .block(i * m_nThrottle,
                   m_nStates * (m_nIter + 1) + m_nJoints * m_ctrlHorizon + i * 2 * m_nThrottle
                       + thrustIdx[0],
                   m_nThrottle,
                   m_nThrottle)
            .setIdentity();
        if (m_counter != (m_ratioSmallLargeStepsPeriod - 1) && i == 0)
        {
            for (int j = 0; j < N_THRUSTS; j++)
            {
                m_lowerBound[j] = (m_robot->getJetThrusts()[j]);
                m_upperBound[j] = (m_robot->getJetThrusts()[j]);
            }
        }
    }
    if (m_counter == (m_ratioSmallLargeStepsPeriod - 1))
    {
        qpInput.setUpdateThrottle(true);
        m_counter = 0;
    } else
    {
        qpInput.setUpdateThrottle(false);
        m_counter++;
    }
    return true;
}

const bool ThrustContraintDD::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
ThrustContraintDD::populateVectorsCollection(QPInput& qpInput,
                                             const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

HankleMatrixConstraint::HankleMatrixConstraint(const int nVar,
                                               const int nStates,
                                               const int horizonLenghtHankleMatrix,
                                               const int nJets,
                                               const int nG,
                                               const int gInitPosition,
                                               const int throttleInitPosition,
                                               const int thrustInitPosition,
                                               const std::vector<std::vector<double>>& inputData,
                                               const std::vector<std::vector<double>>& outputData)
    : IQPConstraint(nVar, 2 * nJets * horizonLenghtHankleMatrix)
{
    m_nStates = nStates;
    m_inputData = inputData;
    m_outputData = outputData;
    m_horizonLenghtHankleMatrix = horizonLenghtHankleMatrix;
    m_nG = nG;
    m_initGPosition = gInitPosition;
    m_initThrottlePosition = throttleInitPosition;
    m_initThrustPosition = thrustInitPosition;
}

const bool HankleMatrixConstraint::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("nIterSmall", m_nsmallSteps))
    {
        yError() << "Parameter 'nIterSmall' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("controlHorizon", m_ctrlHorizonInput))
    {
        yError() << "Parameter 'controlHorizon' not found in the config file.";
        return false;
    }
    double periodMPCLargeSteps;
    if (!ptr->getParameter("periodMPCLargeSteps", periodMPCLargeSteps))
    {
        yError() << "Parameter 'periodMPCLargeSteps' not found in the config file.";
        return false;
    }
    double periodMPCSmallSteps;
    if (!ptr->getParameter("periodMPCSmallSteps", periodMPCSmallSteps))
    {
        yError() << "Parameter 'periodMPCSmallSteps' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("updateHankelMatrix", m_updateHankleMatrix))
    {
        yError() << "Parameter 'updateHankelMatrix' not found in the config file.";
        return false;
    }
    m_ratioSmallLargeStepsPeriod = round(periodMPCLargeSteps / periodMPCSmallSteps);
    m_initDataLength = m_horizonLenghtHankleMatrix - (m_ctrlHorizonInput - m_nsmallSteps + 1);
    m_robot = qpInput.getRobot();
    m_jetModel = qpInput.getJetModel();
    m_nJets = 4;
    return true;
}

void HankleMatrixConstraint::configureDynVectorsSize(QPInput& qpInput)
{
    m_hankleMatrixInput.resize(m_robot->getNJets());
    m_hankleMatrixOutput.resize(m_robot->getNJets());
    m_numCols = m_inputData[0].size() - m_horizonLenghtHankleMatrix;
    for (int i = 0; i < m_robot->getNJets(); i++)
    {
        m_hankleMatrixInput[i].resize(m_horizonLenghtHankleMatrix, m_numCols);
        m_hankleMatrixOutput[i].resize(m_horizonLenghtHankleMatrix, m_numCols);
    }
    this->updateHankleMatrix();
    m_linearMatrix.setZero();
    m_lowerBound.setZero();
    m_upperBound.setZero();
    m_initialThrottleValuesVector.resize(m_nJets);
    m_initialThrustValuesVector.resize(m_nJets);
    for (int i = 0; i < m_nJets; i++)
    {
        m_initialThrottleValuesVector[i].resize(m_initDataLength);
        m_initialThrustValuesVector[i].resize(m_initDataLength);
        for (int j = 0; j < m_initDataLength; j++)
        {
            m_initialThrottleValuesVector[i][j] = 0.0;
            m_initialThrustValuesVector[i][j] = 0.0;
        }
    }
    m_counter = m_ratioSmallLargeStepsPeriod - 1;
}

const bool HankleMatrixConstraint::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    if (m_firstIteriation)
    {
        for (int i = 0; i < m_nJets; i++)
        {
            m_linearMatrix.block(i * m_horizonLenghtHankleMatrix,
                                 m_initGPosition + i * m_nG,
                                 m_horizonLenghtHankleMatrix,
                                 m_numCols)
                = m_hankleMatrixInput[i];
            m_linearMatrix.block(m_nJets * m_horizonLenghtHankleMatrix
                                     + i * m_horizonLenghtHankleMatrix,
                                 m_initGPosition + i * m_nG,
                                 m_horizonLenghtHankleMatrix,
                                 m_numCols)
                = m_hankleMatrixOutput[i];
            // set the constraints for the slack variables
            m_linearMatrix
                .block(m_nJets * m_horizonLenghtHankleMatrix + i * m_horizonLenghtHankleMatrix,
                       m_initGPosition + m_nJets * m_nG + i * m_initDataLength,
                       m_initDataLength,
                       m_initDataLength)
                .setIdentity();
        }
        for (int i = 0; i < m_nJets; i++)
        {
            for (int j = 0; j < (m_ctrlHorizonInput - m_nsmallSteps + 1); j++)
            {
                m_linearMatrix(i * m_horizonLenghtHankleMatrix + m_initDataLength + j,
                               m_initThrottlePosition + j * m_nJets + i)
                    = -1;
                m_linearMatrix(m_nJets * m_horizonLenghtHankleMatrix
                                   + i * m_horizonLenghtHankleMatrix + m_initDataLength + j,
                               m_initThrustPosition + j * m_nJets + i)
                    = -1;
            }
        }
        m_firstIteriation = false;
    }
    if (m_counter == (m_ratioSmallLargeStepsPeriod - 1))
    {
        for (int i = 0; i < m_nJets; i++)
        {
            // update the constraints related to the hankel matrix
            m_linearMatrix.block(i * m_horizonLenghtHankleMatrix,
                                 m_initGPosition + i * m_nG,
                                 m_horizonLenghtHankleMatrix,
                                 m_numCols)
                = m_hankleMatrixInput[i];
            m_linearMatrix.block(m_nJets * m_horizonLenghtHankleMatrix
                                     + i * m_horizonLenghtHankleMatrix,
                                 m_initGPosition + i * m_nG,
                                 m_horizonLenghtHankleMatrix,
                                 m_numCols)
                = m_hankleMatrixOutput[i];
            // remove the first element of the vectors
            m_initialThrottleValuesVector[i].erase(m_initialThrottleValuesVector[i].begin());
            m_initialThrustValuesVector[i].erase(m_initialThrustValuesVector[i].begin());
            // add the new element at the end of the vectors
            m_initialThrottleValuesVector[i].push_back(
                m_jetModel->compute_u2v(qpInput.getThrottleMPC()[i]));
            m_initialThrustValuesVector[i].push_back(m_robot->getJetThrusts()[i]
                                                     - qpInput.getOutputQPThrust()[i]);
        }
        for (int i = 0; i < m_nJets; i++)
        {
            for (int j = 0; j < m_initDataLength; j++)
            {
                m_upperBound[i * m_horizonLenghtHankleMatrix + j]
                    = m_initialThrottleValuesVector[i][j];
                m_lowerBound[i * m_horizonLenghtHankleMatrix + j]
                    = m_initialThrottleValuesVector[i][j];
                m_upperBound[m_nJets * m_horizonLenghtHankleMatrix + i * m_horizonLenghtHankleMatrix
                             + j]
                    = m_initialThrustValuesVector[i][j];
                m_lowerBound[m_nJets * m_horizonLenghtHankleMatrix + i * m_horizonLenghtHankleMatrix
                             + j]
                    = m_initialThrustValuesVector[i][j];
            }
        }
        m_counter = 0;
    } else
    {
        m_counter++;
    }

    return true;
}

void HankleMatrixConstraint::updateHankleMatrix()
{
    for (int i = 0; i < m_inputData.size(); i++)
    {
        for (int j = 0; j < m_horizonLenghtHankleMatrix; j++)
        {
            for (int k = 0; k < m_numCols; k++)
            {
                m_hankleMatrixInput[i](j, k) = m_jetModel->compute_u2v(m_inputData[i][k + j]);
                m_hankleMatrixOutput[i](j, k) = (m_outputData[i][k + j]);
            }
        }
    }
    // Check if the Hankle matrix is persistently exciting
    for (int i = 0; i < m_nJets; i++)
    {
        if (!isHankleMatrixPersistelyExciting(m_hankleMatrixInput[i]))
        {
            yError() << "Input Hankle matrix is not persistently exciting for jet " << i
                     << ". Please check the input data.";
        }
        if (!isHankleMatrixPersistelyExciting(m_hankleMatrixOutput[i]))
        {
            yError() << "Output Hankle matrix is not persistently exciting for jet " << i
                     << ". Please check the output data.";
        }
    }
}

bool HankleMatrixConstraint::isHankleMatrixPersistelyExciting(Eigen::MatrixXd& hankelMatrix,
                                                              const double tolerance)
{
    int max_possible_rank = std::min(hankelMatrix.rows(), hankelMatrix.cols());
    int order = 10;

    if (order > hankelMatrix.rows())
    {
        std::cerr << "Warning: Requested order (" << order
                  << ") for PE check exceeds Hankel matrix row dimension (" << hankelMatrix.rows()
                  << "). This check might not be meaningful." << std::endl;
        // Depending on interpretation, you might throw or just return false here
    }

    if (tolerance <= 0.0)
    {
        throw std::invalid_argument("Tolerance for PE check must be positive.");
    }

    // Handle empty matrix case
    if (hankelMatrix.rows() == 0 || hankelMatrix.cols() == 0)
    {
        // An empty matrix has rank 0. It's PE only if the required rank is 0.
        return (order == 0);
    }

    // --- Use ColPivHouseholderQR for Rank Computation ---
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(hankelMatrix);

    // Get the absolute values of the diagonal elements of the R factor.
    // The rank is the number of diagonal elements whose magnitude is > tolerance.
    // Eigen's internal representation might store R differently, but diagonal() on matrixQR() gives
    // access. Size of diagonal is min(rows, cols).
    const auto& R_diag_abs = qr.matrixQR().diagonal().cwiseAbs();

    // Count the number of diagonal elements above the *absolute* tolerance

    int computed_rank = 0;
    for (int i = 0; i < R_diag_abs.size(); ++i)
    {
        if (R_diag_abs(i) > tolerance)
        {
            computed_rank++;
        }
    }

    bool isPE = (computed_rank >= order);
    if (!isPE)
    {
        yError() << "Hankel matrix is not persistently exciting. Computed rank: " << computed_rank
                 << ", required rank: " << order;
    }

    // Check if the computed rank meets the required rank 'order'
    return isPE;
}

const bool HankleMatrixConstraint::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
HankleMatrixConstraint::populateVectorsCollection(QPInput& qpInput,
                                                  const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

ArtificialEquilibriumConstraint::ArtificialEquilibriumConstraint(
    const int nVar,
    const int nStates,
    const int nIter,
    const int nArtificialEquilibriumStates,
    const int ArtificialEquilibriumStatesInitPosition)
    : IQPConstraint(nVar, nArtificialEquilibriumStates)
{
    m_ArtificialEquilibriumStatesInitPosition = ArtificialEquilibriumStatesInitPosition;
    m_nArtificialEquilibriumStates = nArtificialEquilibriumStates;
    m_nIter = nIter;
    m_nStates = nStates;
    std::cout << "ArtificialEquilibriumConstraint created with nArtificialEquilibriumStates: "
              << m_nArtificialEquilibriumStates << std::endl;
}

const bool ArtificialEquilibriumConstraint::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    return true;
}

void ArtificialEquilibriumConstraint::configureDynVectorsSize(QPInput& qpInput)
{
}

const bool ArtificialEquilibriumConstraint::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    if (m_firstIteration)
    {
        int constraintIdx = 0;
        // CoM_{L} = x_{s, com}
        m_linearMatrix
            .block(constraintIdx,
                   m_nStates * m_nIter + CoMPosIdx[0],
                   CoMPosIdx.size(),
                   CoMPosIdx.size())
            .setIdentity();
        m_linearMatrix.block(constraintIdx,
                             m_ArtificialEquilibriumStatesInitPosition,
                             CoMPosIdx.size(),
                             CoMPosIdx.size())
            = -Eigen::MatrixXd::Identity(CoMPosIdx.size(), CoMPosIdx.size());
        constraintIdx += CoMPosIdx.size();
        // linMom_{L} = x_{s, linMom}
        m_linearMatrix
            .block(constraintIdx,
                   m_nStates * m_nIter + linMomIdx[0],
                   linMomIdx.size(),
                   linMomIdx.size())
            .setIdentity();
        m_linearMatrix.block(constraintIdx,
                             m_ArtificialEquilibriumStatesInitPosition + constraintIdx,
                             linMomIdx.size(),
                             linMomIdx.size())
            = -Eigen::MatrixXd::Identity(linMomIdx.size(), linMomIdx.size());
        constraintIdx += linMomIdx.size();
        // rpy_{L} = x_{s, rpy}
        m_linearMatrix
            .block(constraintIdx, m_nStates * m_nIter + rpyIdx[0], rpyIdx.size(), rpyIdx.size())
            .setIdentity();
        m_linearMatrix.block(constraintIdx,
                             m_ArtificialEquilibriumStatesInitPosition + constraintIdx,
                             rpyIdx.size(),
                             rpyIdx.size())
            = -Eigen::MatrixXd::Identity(rpyIdx.size(), rpyIdx.size());
        constraintIdx += rpyIdx.size();
        // angMom_{L} = x_{s, angMom}
        m_linearMatrix
            .block(constraintIdx,
                   m_nStates * m_nIter + angMomIdx[0],
                   angMomIdx.size(),
                   angMomIdx.size())
            .setIdentity();
        m_linearMatrix.block(constraintIdx,
                             m_ArtificialEquilibriumStatesInitPosition + constraintIdx,
                             angMomIdx.size(),
                             angMomIdx.size())
            = -Eigen::MatrixXd::Identity(angMomIdx.size(), angMomIdx.size());
        m_lowerBound.setZero();
        m_upperBound.setZero();

        m_firstIteration = false;
    }
    return true;
}

const bool ArtificialEquilibriumConstraint::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool ArtificialEquilibriumConstraint::populateVectorsCollection(
    QPInput& qpInput, const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

ThrustHatConstraint::ThrustHatConstraint(const int nVar,
                                         const double periodMPC,
                                         const int nThrustHat,
                                         const int thrustHatInitPosition,
                                         const int thrustHatDotInitPosition)
    : IQPConstraint(nVar, 2 * 4)
{
    m_nThrustHat = nThrustHat;
    m_thrustHatInitPosition = thrustHatInitPosition;
    m_periodMPC = periodMPC;
    m_thrustHatDotInitPosition = thrustHatDotInitPosition;
}

const bool ThrustHatConstraint::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    m_robot = qpInput.getRobot();
    m_jetModel = qpInput.getJetModel();
    return true;
}

void ThrustHatConstraint::configureDynVectorsSize(QPInput& qpInput)
{
    m_nJets = m_robot->getNJets();
}

const bool ThrustHatConstraint::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    m_linearMatrix.setZero();
    m_lowerBound.setZero();
    m_upperBound.setZero();
    for (int j = 0; j < m_robot->getNJets(); j++)
    {
        m_thrust = m_robot->getJetThrusts()[j];
        m_thrustDot = qpInput.getEstimatedThrustDot()[j];
        m_linearMatrix(j, m_thrustHatInitPosition + j + thrustHatIdx[0]) = 1;
        m_linearMatrix(j + m_robot->getNJets(), m_thrustHatInitPosition + j + thrustHatDotIdx[0])
            = 1;
        m_lowerBound(j) = m_thrust;
        m_upperBound(j) = m_thrust;
        m_lowerBound(j + m_robot->getNJets()) = m_thrustDot;
        m_upperBound(j + m_robot->getNJets()) = m_thrustDot;
    }

    return true;
}

const bool ThrustHatConstraint::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
ThrustHatConstraint::populateVectorsCollection(QPInput& qpInput,
                                               const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

ThrustSummationConstraint::ThrustSummationConstraint(const int nVar,
                                                     const int nStates,
                                                     const int nIter,
                                                     const int nSmallSteps,
                                                     const int ctrlHorizon,
                                                     const int thrustInitPosition,
                                                     const int thrustHatInitPosition,
                                                     const int thrustTildeInitPosition)
    : IQPConstraint(nVar, N_THRUSTS * (ctrlHorizon - nSmallSteps + 1))
{
    m_thrustTildeInitPosition = thrustTildeInitPosition;
    m_thrustHatInitPosition = thrustHatInitPosition;
    m_thrustInitPosition = thrustInitPosition;
    m_nIter = nIter;
    m_nStates = nStates;
    m_nSmallSteps = nSmallSteps;
    m_ctrlHorizon = ctrlHorizon;
}

const bool ThrustSummationConstraint::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    return true;
}

void ThrustSummationConstraint::configureDynVectorsSize(QPInput& qpInput)
{
    m_robot = qpInput.getRobot();
    m_nJets = m_robot->getNJets();
    m_nJoints = deltaJointIdx.size();
    m_nThrottle = throttleIdx.size();
}

const bool ThrustSummationConstraint::computeConstraintsMatrixAndBounds(QPInput& qpInput)
{
    m_linearMatrix.setZero();
    for (int i = 0; i < m_ctrlHorizon - m_nSmallSteps + 1; i++)
    {
        for (int j = 0; j < m_nJets; j++)
        {
            m_linearMatrix(i * m_nJets + j, m_thrustInitPosition + i * m_nJets + j) = 1;
            m_linearMatrix(i * m_nJets + j,
                           m_thrustHatInitPosition + (i + 1) * 2 * m_nJets + j + thrustHatIdx[0])
                = -1;
            m_linearMatrix(i * m_nJets + j, m_thrustTildeInitPosition + i * m_nJets + j) = -1;
        }
    }
    m_lowerBound.setZero();
    m_upperBound.setZero();
    return true;
}

const bool ThrustSummationConstraint::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
ThrustSummationConstraint::populateVectorsCollection(QPInput& qpInput,
                                                     const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

JetDynamicsConstraint::JetDynamicsConstraint(const unsigned int nVar,
                                             const unsigned int nStates,
                                             const unsigned int nInput,
                                             const unsigned int nIter,
                                             const unsigned int startIdx)
    : IQPConstraintMPCDynamic(nVar, nStates, nInput, nIter, startIdx)
{
    m_nIter = nIter;
    m_nStates = nStates;
    m_jetDynamic = std::make_shared<JetDynamicVS>(nStates, nInput);
}

const bool JetDynamicsConstraint::readConfigParameters(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    auto ptr = parametersHandler.lock();
    m_jetDynamic->configure(ptr, qpInput);
    return true;
}

void JetDynamicsConstraint::configureDynVectorsSize(QPInput& qpInput)
{
}

bool JetDynamicsConstraint::updateDynamicConstraints(QPInput& qpInput)
{
    m_jetDynamic->updateInitialStates(qpInput);
    m_jetDynamic->getAMatrixDynamicClass(m_A);
    m_jetDynamic->getBMatrixDynamicClass(m_B);
    m_jetDynamic->getCVectorDynamicClass(m_c);
    return true;
}

const bool JetDynamicsConstraint::configureVectorsCollectionServer(QPInput& qpInput)
{
    return true;
}

const bool
JetDynamicsConstraint::populateVectorsCollection(QPInput& qpInput,
                                                 const Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    return true;
}

} // namespace DDMPC
