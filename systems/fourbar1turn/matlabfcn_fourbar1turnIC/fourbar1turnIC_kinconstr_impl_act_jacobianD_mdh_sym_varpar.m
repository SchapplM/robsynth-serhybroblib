% Jacobian time derivative of explicit kinematic constraints of
% fourbar1turnIC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = fourbar1turnIC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:11
% EndTime: 2020-05-07 11:33:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (6->4), mult. (8->6), div. (0->0), fcn. (4->4), ass. (0->4)
t6 = pkin(3) * (qJD(2) + qJD(3));
t5 = pkin(2) * qJD(2);
t4 = qJ(2) + qJ(3);
t1 = [0, -cos(t4) * t6 + cos(qJ(2)) * t5; 0, -sin(t4) * t6 + sin(qJ(2)) * t5; 0, 0;];
PhiD_a = t1;
