#include <FlightControlUtils.h>
#include <dataDrivenMPC/DDconstant.h>
#include <dataDrivenMPC/constraintsDDMPC.h>
#include <dataDrivenMPC/costsDDMPC.h>
#include <dataDrivenMPC/dataDrivenMPC.h>

const bool DataDrivenMPC::setCostAndConstraints(
    std::weak_ptr<BipedalLocomotion::ParametersHandler::IParametersHandler> parametersHandler,
    QPInput& qpInput)
{
    if (!m_hankleMatrixSet)
    {
        yError() << "DataFusedMPC::setCostAndConstraints: the Hankel matrices are "
                    "not set; please first call setHankleMatrices()";
        return false;
    }
    auto ptr = parametersHandler.lock();
    if (!ptr->getParameter("controlledJoints", m_controlledJoints))
    {
        yError() << "Parameter 'controlledJoints' not found in the config file.";
        return false;
    }
    m_nCtrlJoints = m_controlledJoints.size();
    if (m_nCtrlJoints != deltaJointIdx.size())
    {
        yError() << "The number of controlled joints defined in the systemDynamic.h file is "
                    "different from the size of the 'controlledJoints' parameter";
        return false;
    }
    if (!ptr->getParameter("nIter", m_nIter))
    {
        yError() << "Parameter 'nIter' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("nIterSmall", m_nIterSmall))
    {
        yError() << "Parameter 'nIterSmall' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("controlHorizon", m_ctrlHorizon))
    {
        yError() << "Parameter 'controlHorizon' not found in the config file.";
        return false;
    }
    if (!ptr->getParameter("HankleMatrixHorizon", m_horizonLenghtHankleMatrix))
    {
        yError() << "Parameter 'HankleMatrixHorizon' not found in the config file.";
        return false;
    }
    bool useResidualThrust = true;
    m_robot = qpInput.getRobot();
    m_nJets = m_robot->getNJets();
    m_jetModel = qpInput.getJetModel();
    m_nStates = angMomIdx[2] + 1;
    m_nInput = m_nCtrlJoints + m_nJets;
    int gParamNumber = m_inputData[0].size() - m_horizonLenghtHankleMatrix;
    // thurst hat is the thrust coming from the second order model
    int nThrustHat = m_nJets * (m_ctrlHorizon - m_nIterSmall + 2);
    int nThrustHatDot = m_nJets * (m_ctrlHorizon - m_nIterSmall + 2);
    // thrust tilde is the thrust residual coming from the data driven model
    int nThrustTilde = m_nJets * (m_ctrlHorizon - m_nIterSmall + 1);

    // number of slack variables
    int slackVarNumber
        = m_nJets * (m_horizonLenghtHankleMatrix - (m_ctrlHorizon - m_nIterSmall + 1));
    if (!useResidualThrust)
    {
        slackVarNumber = 0;
        gParamNumber = 0;
    }

    // number of artificial equilibrium states variables
    m_nArtificialEquilibriumStates = 0;

    // optimization variable vector:
    // z = [(x, h^l, phi, h^w)_i, delta_s, T, g, slack, artificialEq, (T_hat, T_hat_dot)_i, v,
    // T_tilde]

    // number of variables
    m_nVar = m_nStates * (m_nIter + 1) + m_nCtrlJoints * m_ctrlHorizon
             + 2 * m_nJets * (m_ctrlHorizon - m_nIterSmall + 1) + m_nJets * gParamNumber
             + slackVarNumber + m_nArtificialEquilibriumStates + nThrustHat + nThrustHatDot
             + nThrustTilde;
    m_gParamInitPosition = m_nStates * (m_nIter + 1) + m_nCtrlJoints * m_ctrlHorizon
                           + m_nJets * (m_ctrlHorizon - m_nIterSmall + 1);
    int thrustInitPosition = m_nStates * (m_nIter + 1) + m_nCtrlJoints * m_ctrlHorizon;
    int slackVarInitPosition = m_nStates * (m_nIter + 1) + m_nCtrlJoints * m_ctrlHorizon
                               + m_nJets * (m_ctrlHorizon - m_nIterSmall + 1)
                               + m_nJets * gParamNumber;
    m_artificialEquilibriumStatesInitPosition
        = m_nStates * (m_nIter + 1) + m_nCtrlJoints * m_ctrlHorizon
          + m_nJets * (m_ctrlHorizon - m_nIterSmall + 1) + m_nJets * gParamNumber + slackVarNumber;
    m_thrustHatInitPosition
        = m_artificialEquilibriumStatesInitPosition + m_nArtificialEquilibriumStates;
    m_thrustHatDotInitPosition = m_thrustHatInitPosition + nThrustHat;
    m_throttleInitPosition = m_thrustHatDotInitPosition + nThrustHatDot;
    int thrustTildeInitPosition = m_throttleInitPosition + nThrustTilde;
    std::cout << "thrust init position: " << thrustInitPosition
              << " throttle init position: " << m_throttleInitPosition << " throttle end position: "
              << m_throttleInitPosition + m_nJets * (m_ctrlHorizon - m_nIterSmall + 1) << std::endl;
    std::cout << " thrust tilde init position: " << thrustTildeInitPosition << std::endl;
    m_jointSelectorVector.clear();
    for (auto joint : m_controlledJoints)
    {
        for (int i = 0; i < m_robot->getNJoints(); i++)
        {
            if (joint == m_robot->getJointName(i))
            {
                m_jointSelectorVector.emplace_back(i);
            }
        }
    }

    // resize the vectors
    m_jointsPositionReference.resize(m_robot->getNJoints());
    m_jointsPositionReference = m_robot->getJointPos();
    m_previousState.resize(m_nStates);
    m_QPSolution.resize(m_nVar);
    m_deltaJointsPositionReference.resize(m_nCtrlJoints);
    m_thrustReference.resize(m_nJets);
    m_throttleReference.resize(m_nJets);
    m_statesSolution.resize(m_nStates * (m_nIter + 1));
    m_inputSolution.resize(m_nCtrlJoints * m_nIter + m_nJets * (m_nIter - m_nIterSmall + 1));

    m_vectorCosts.emplace_back(
        std::make_unique<DDMPC::ReferenceTrackingCost>(m_nVar, m_nStates, m_nIter));
    m_vectorCosts.emplace_back(
        std::make_unique<DDMPC::RegualarizationCost>(m_nVar,
                                                     m_nStates,
                                                     m_nCtrlJoints,
                                                     m_nJets,
                                                     m_throttleInitPosition));
    m_vectorCosts.emplace_back(
        std::make_unique<DDMPC::ThrottleInitialValueCost>(m_nVar,
                                                          m_nStates,
                                                          m_nCtrlJoints,
                                                          m_nJets,
                                                          m_throttleInitPosition));
    m_vectorCosts.emplace_back(
        std::make_unique<DDMPC::RegularizationGParametersCost>(m_nVar,
                                                               m_gParamInitPosition,
                                                               m_nJets * gParamNumber));
    m_vectorCosts.emplace_back(
        std::make_unique<DDMPC::RegularizationSlackVariableCost>(m_nVar,
                                                                 slackVarInitPosition,
                                                                 slackVarNumber));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::ConstraintSystemDynamicDD>(m_nVar,
                                                           m_nStates,
                                                           m_nCtrlJoints,
                                                           m_nJets,
                                                           m_nIter));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::ConstraintInitialStateDD>(m_nStates, m_nVar));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::ThrottleConstraintDD>(m_nVar,
                                                      m_nStates,
                                                      m_nIter,
                                                      m_nIterSmall,
                                                      m_ctrlHorizon,
                                                      m_throttleInitPosition));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::JointPositionConstraintDD>(m_nVar,
                                                           m_nStates,
                                                           m_nCtrlJoints,
                                                           m_nIter,
                                                           m_ctrlHorizon));
    // m_vectorConstraints.emplace_back(
    //     std::make_unique<DDMPC::ThrustContraintDD>(m_nVar,
    //                                                m_nStates,
    //                                                m_nIter,
    //                                                m_nIterSmall,
    //                                                m_ctrlHorizon,
    //                                                thrustHatInitPosition));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::HankleMatrixConstraint>(m_nVar,
                                                        m_nStates,
                                                        m_horizonLenghtHankleMatrix,
                                                        m_nJets,
                                                        gParamNumber,
                                                        m_gParamInitPosition,
                                                        m_throttleInitPosition,
                                                        thrustTildeInitPosition,
                                                        m_inputData,
                                                        m_outputData));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::ThrustHatConstraint>(m_nVar,
                                                     0.1,
                                                     nThrustHat,
                                                     m_thrustHatInitPosition,
                                                     m_throttleInitPosition));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::ThrustSummationConstraint>(m_nVar,
                                                           m_nStates,
                                                           m_nIter,
                                                           m_nIterSmall,
                                                           m_ctrlHorizon,
                                                           thrustInitPosition,
                                                           m_thrustHatInitPosition,
                                                           thrustTildeInitPosition));
    m_vectorConstraints.emplace_back(
        std::make_unique<DDMPC::JetDynamicsConstraint>(m_nVar,
                                                       2 * m_robot->getNJets(),
                                                       m_robot->getNJets(),
                                                       m_ctrlHorizon - m_nIterSmall + 1,
                                                       m_thrustHatInitPosition));
    return true;
}

const bool DataDrivenMPC::solveMPC()
{
    this->solve();
    if (this->getQPProblemStatus() == OsqpEigen::Status::Solved)
    {
        m_QPSolution = this->getSolution();

        // extract the solution
        m_statesSolution = m_QPSolution.head(m_nStates * (m_nIter + 1));
        m_inputSolution = m_QPSolution.segment(m_nStates * (m_nIter + 1),
                                               m_nCtrlJoints * m_ctrlHorizon
                                                   + m_nJets * (m_ctrlHorizon - m_nIterSmall + 1));
        m_deltaJointsPositionReference = m_inputSolution.segment(0, m_nCtrlJoints);
        m_throttleReference = m_QPSolution.segment(m_throttleInitPosition, m_nJets);
        m_thrustReference = m_inputSolution.segment(m_nCtrlJoints * m_ctrlHorizon, m_nJets);
        m_finalState = m_statesSolution.tail(m_nStates);
        for (int i = 0; i < m_jointSelectorVector.size(); i++)
        {
            m_jointsPositionReference(m_jointSelectorVector[i])
                += m_deltaJointsPositionReference(i);
        }
        auto lastThrustHat = m_QPSolution.segment(m_thrustHatDotInitPosition - m_robot->getNJets(),
                                                  thrustHatIdx.size());
    }

    return true;
}

const double DataDrivenMPC::getValueFunction()
{
    if (this->getQPProblemStatus() == OsqpEigen::Status::Solved)
    {
        return this->getCostValue();
    } else
    {
        yError() << "DataFusedMPC::getValueFunction: the MPC problem is not "
                    "solved";
        return 0.0;
    }
}

const bool DataDrivenMPC::setHankleMatrices(const std::vector<std::vector<double>>& inputData,
                                            const std::vector<std::vector<double>>& outputData)
{
    m_hankleMatrixSet = true;
    if ((inputData.size() != N_THRUSTS) || (outputData.size() != N_THRUSTS))
    {
        yError() << "DataFusedMPC::setHankleMatrices: the input and output vector "
                    "must have the size: "
                 << N_THRUSTS << " instead of " << inputData.size() << " and " << outputData.size();
        return false;
    }
    m_inputData = inputData;
    m_outputData = outputData;
    return true;
}

const bool DataDrivenMPC::getMPCSolution(Eigen::Ref<Eigen::VectorXd> qpSolution)
{
    if (qpSolution.size() != m_nInput)
    {
        yError() << "DataFusedMPC::getMPCSolution: wrong size of the input vector";
        return false;
    }
    qpSolution = m_inputSolution;
    return true;
}

const bool DataDrivenMPC::getJointsReferencePosition(Eigen::Ref<Eigen::VectorXd> jointsPosition)
{
    if (jointsPosition.size() != m_robot->getNJoints())
    {
        yError() << "DataFusedMPC::getJointsReferencePosition: wrong size of the "
                    "input "
                    "vector";
        return false;
    }
    jointsPosition = m_jointsPositionReference;
    return true;
}

const bool DataDrivenMPC::getThrottleReference(Eigen::Ref<Eigen::VectorXd> throttle)
{
    if (throttle.size() != m_nJets)
    {
        yError() << "DataFusedMPC::getThrottleReference: wrong size of the input "
                    "vector";
        return false;
    }

    for (int i = 0; i < m_nJets; i++)
    {
        throttle(i) = m_jetModel->destandardizeThrottle_u2T(m_throttleReference(i));
    }
    return true;
}

const bool DataDrivenMPC::getThrustReference(Eigen::Ref<Eigen::VectorXd> thrust)
{
    if (thrust.size() != m_nJets)
    {
        yError() << "DataFusedMPC::getThrustReference: wrong size of the input "
                    "vector";
        return false;
    }
    for (int i = 0; i < m_nJets; i++)
    {
        thrust(i) = (m_thrustReference(i));
    }
    return true;
}

const bool DataDrivenMPC::getFinalCoMPosition(Eigen::Ref<Eigen::VectorXd> finalCoMPosition)
{
    if (finalCoMPosition.size() != 3)
    {
        yError() << "DataFusedMPC::getFinalCoMPosition: wrong size of the input "
                    "vector";
        return false;
    }
    finalCoMPosition = m_finalState.segment(CoMPosIdx[0], CoMPosIdx.size());
    return true;
}

const bool DataDrivenMPC::getFinalLinMom(Eigen::Ref<Eigen::VectorXd> finalLinMom)
{
    if (finalLinMom.size() != 3)
    {
        yError() << "DataFusedMPC::getFinalLinMom: wrong size of the input vector";
        return false;
    }
    finalLinMom = m_finalState.segment(linMomIdx[0], linMomIdx.size());
    return true;
}

const bool DataDrivenMPC::getFinalRPY(Eigen::Ref<Eigen::VectorXd> finalRPY)
{
    if (finalRPY.size() != 3)
    {
        yError() << "DataFusedMPC::getFinalRPY: wrong size of the input vector";
        return false;
    }
    finalRPY = m_finalState.segment(rpyIdx[0], rpyIdx.size());
    return true;
}

const bool DataDrivenMPC::getFinalAngMom(Eigen::Ref<Eigen::VectorXd> finalAngMom)
{
    if (finalAngMom.size() != 3)
    {
        yError() << "DataFusedMPC::getFinalAngMom: wrong size of the input vector";
        return false;
    }
    finalAngMom = m_finalState.segment(angMomIdx[0], angMomIdx.size());
    return true;
}

const bool DataDrivenMPC::getThrustHat(Eigen::Ref<Eigen::VectorXd> thrustHat)
{
    if (thrustHat.size() != m_robot->getNJets())
    {
        yError() << "DataFusedMPC::getThrustHat: wrong size of the input vector";
        return false;
    }
    thrustHat = m_QPSolution.segment(m_thrustHatInitPosition + 2 * m_robot->getNJets(),
                                     thrustHatIdx.size());
    return true;
}

double DataDrivenMPC::getNStatesMPC() const
{
    return m_nStates;
}

double DataDrivenMPC::getNInputMPC() const
{
    return m_nInput;
}
