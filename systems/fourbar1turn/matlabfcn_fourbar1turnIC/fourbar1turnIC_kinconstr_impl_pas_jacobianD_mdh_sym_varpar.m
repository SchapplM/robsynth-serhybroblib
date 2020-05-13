% Jacobian time derivative of explicit kinematic constraints of
% fourbar1turnIC
% with respect to passive joint coordinates
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
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = fourbar1turnIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:12
% EndTime: 2020-05-07 11:33:12
% DurationCPUTime: 0.01s
% Computational Cost: add. (8->6), mult. (8->6), div. (0->0), fcn. (4->4), ass. (0->4)
t11 = pkin(3) * (qJD(2) + qJD(3));
t10 = pkin(4) * qJD(4);
t9 = qJ(2) + qJ(3);
t1 = [-cos(t9) * t11, -cos(qJ(4)) * t10, 0; -sin(t9) * t11, -sin(qJ(4)) * t10, 0; 0, 0, 0;];
PhiD_p = t1;
