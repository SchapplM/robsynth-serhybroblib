% Jacobian time derivative of explicit kinematic constraints of
% picker2Dm1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:55
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = picker2Dm1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:54:42
% EndTime: 2020-05-11 05:54:43
% DurationCPUTime: 0.03s
% Computational Cost: add. (23->13), mult. (24->15), div. (0->0), fcn. (12->10), ass. (0->12)
t19 = qJD(1) + qJD(2);
t30 = pkin(3) * t19;
t29 = pkin(4) * (qJD(4) + t19);
t28 = pkin(5) * (qJD(1) + qJD(8));
t27 = pkin(1) * qJD(1);
t26 = pkin(3) * qJD(7);
t21 = qJ(1) + qJ(2);
t25 = sin(qJ(1)) * t27;
t24 = cos(qJ(1)) * t27;
t20 = qJ(1) + qJ(8);
t17 = qJ(4) + t21;
t1 = [-cos(t17) * t29 - cos(t21) * t30 - t24, -sin(qJ(7)) * t26; -sin(t17) * t29 - sin(t21) * t30 - t25, cos(qJ(7)) * t26; -cos(t20) * t28 + t24, 0; -sin(t20) * t28 + t25, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
PhiD_a = t1;
